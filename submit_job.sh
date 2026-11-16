#!/bin/bash --login
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1    
#SBATCH --cpus-per-task=20         
#SBATCH --mem=20G 
#SBATCH --job-name=zinc_docking 
#SBATCH --time=12:00:00 
#SBATCH --qos=normal
#SBATCH --partition=<>
#SBATCH --account=<your_account>
#SBATCH -o slurm-%j.output
#SBATCH -e slurm-%j.error

# =========================================================
# CONFIGURATION
# =========================================================

VINA_SCRIPTS="${VINA_SCRIPTS:-/your/directory/here/docking/scripts/AutoDock-Vina/example/autodock_scripts}"

# Set FORCE=1 (e.g. FORCE=1 sbatch zinc_docking.sh) to redo every step even if
# its output already exists.
FORCE="${FORCE:-0}"

# =========================================================
# STEP 1: SETUP ENVIRONMENT
# =========================================================

source "$(conda info --base)/etc/profile.d/conda.sh"

conda env list | grep -qE '^docking\s' || conda env create -f environment.yml
conda activate docking
set -uo pipefail
shopt -s nullglob       

JOB_DIR="$(pwd)"
[ -f "$VINA_SCRIPTS/zinc_pseudo.py" ] || { echo "ERROR: $VINA_SCRIPTS/zinc_pseudo.py not found. Set VINA_SCRIPTS."; exit 1; }
VINA_SCRIPTS="$(realpath "$VINA_SCRIPTS")"
[ -f AD4Zn.dat ] || { echo "ERROR: AD4Zn.dat not found"; exit 1; }
AD4Zn_src="$(realpath AD4Zn.dat)"

python "$VINA_SCRIPTS/zinc_pseudo.py" -h > /dev/null 2>&1 \
    || { echo "ERROR: zinc_pseudo.py failed to run under $(python --version 2>&1)"; exit 1; }

NCPU="${SLURM_CPUS_PER_TASK:-8}"
FAILED=()                  # collects steps that failed

# up_to_date OUTPUT INPUT...
# True if OUTPUT exists, is non-empty, and is newer than every INPUT (and FORCE=0).
up_to_date() {
    local out="$1"; shift
    { [ "$FORCE" = 0 ] && [ -s "$out" ]; } || return 1
    local f
    for f in "$@"; do
        [ "$out" -nt "$f" ] || return 1
    done
}

# =========================================================
# STEP 2: PREPARE LIGANDS  (skipped per-ligand if already done)
# =========================================================
# Supported formats: PDB, MOL2, SDF
# (PDB lacks bond orders, so SDF/MOL2 inputs are safer.)
declare -A seen_names

for lig in ligands/*.{pdb,mol2,sdf}; do
    base=$(basename "$lig")
    name="${base%.*}"
    lig_dir="ligands/${name}"
    sdf_out="${lig_dir}/${name}_H.sdf"

    # Guard against ligand1.pdb and ligand1.sdf both mapping to ligands/ligand1/
    if [[ -n "${seen_names[$name]:-}" ]]; then
        echo "WARNING: '$base' has the same base name as '${seen_names[$name]}'; skipping it"
        continue
    fi
    seen_names[$name]="$base"

    mkdir -p "$lig_dir"

    # 2a: Add hydrogens / standardize (scrub can emit several protomers/tautomers)
    if up_to_date "$sdf_out" "$lig"; then
        echo "Skipping scrub for $name (already done)"
    else
        echo "Preparing ligand $name..."
        rm -f "$sdf_out" "$lig_dir"/*.pdbqt      # clear stale outputs
        scrub.py "$lig" -o "$sdf_out" && [ -s "$sdf_out" ] \
            || { echo "scrub.py produced no output for $name, skipping"; rm -f "$sdf_out"; FAILED+=("ligand:$name"); continue; }
    fi

    # 2b: Convert to PDBQT (skip if any PDBQT in the folder is newer than the SDF)
    if [ "$FORCE" = 0 ] && [ -n "$(find "$lig_dir" -maxdepth 1 -name '*.pdbqt' -newer "$sdf_out" -print -quit)" ]; then
        echo "Skipping Meeko for $name (already done)"
    else
        # If scrub produced multiple molecules, write one PDBQT per molecule
        # (file names then come from the molecule titles).
        n_mols=$(grep -c '^\$\$\$\$' "$sdf_out")
        if [ "$n_mols" -gt 1 ]; then
            echo "  $name produced $n_mols molecules; writing one PDBQT each"
            mk_prepare_ligand.py -i "$sdf_out" --multimol_outdir "$lig_dir"
        else
            mk_prepare_ligand.py -i "$sdf_out" -o "${lig_dir}/${name}.pdbqt"
        fi || { echo "mk_prepare_ligand.py failed for $name, skipping"; FAILED+=("ligand:$name"); continue; }
    fi
done

# make_gpf RECEPTOR_PDBQT LIGAND_PDBQT OUT_GPF CENTER NPTS
# Python 3 replacement for prepare_gpf4zn.py (which needs MGLTools/Python 2).
# Writes an AutoGrid4 GPF with the AutoDock4Zn potentials; receptor and ligand
# atom types are read straight from the two PDBQT files.
make_gpf() {
    python - "$1" "$2" "$3" "$4" "$5" 0.375 AD4Zn.dat << 'PY'
import sys, os

rec, lig, out, center, npts, spacing, param = sys.argv[1:8]

def atom_types(path):
    seen = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                t = line.split()[-1]
                if t not in seen:
                    seen.append(t)
    return seen

prefix = os.path.basename(rec)
if prefix.endswith(".pdbqt"):
    prefix = prefix[:-6]
rtypes, ltypes = atom_types(rec), atom_types(lig)
cx, cy, cz = center.split(",")
nx, ny, nz = npts.split(",")

lines = [
    f"npts {nx} {ny} {nz}",
    f"parameter_file {param}",
    f"gridfld {prefix}.maps.fld",
    f"spacing {spacing}",
    f"receptor_types {' '.join(rtypes)}",
    f"ligand_types {' '.join(ltypes)}",
    f"receptor {os.path.basename(rec)}",
    f"gridcenter {cx} {cy} {cz}",
    "smooth 0.5",
]
lines += [f"map {prefix}.{t}.map" for t in ltypes]
lines += [
    f"elecmap {prefix}.e.map",
    f"dsolvmap {prefix}.d.map",
    "dielectric -0.1465",
    # AutoDock4Zn pairwise potentials (same as prepare_gpf4zn.py)
    "nbp_r_eps 0.25 23.2135 12 6 NA TZ",
    "nbp_r_eps 2.1   3.8453 12 6 OA Zn",
    "nbp_r_eps 2.25  7.5914 12 6 SA Zn",
    "nbp_r_eps 1.0   0.0    12 6 HD Zn",
    "nbp_r_eps 2.0   0.0060 12 6 NA Zn",
    "nbp_r_eps 2.0   0.2966 12 6  N Zn",
]
with open(out, "w") as fh:
    fh.write("\n".join(lines) + "\n")
PY
}

# =========================================================
# Helper: grid + dock one receptor/ligand pair (runs in a subshell)
# =========================================================
dock_pair() {
    local name="$1" lig_abs="$2" lig_base="$3" tz_src="$4" pair_dir="$5" center="$6"
    local tz_file="${name}_tz.pdbqt"
    local gpf_file="${name}_tz+${lig_base}.gpf"
    local glg_file="${name}_tz+${lig_base}.glg"
    local out_file="${name}_tz+${lig_base}.pdbqt"

    cp -f "$AD4Zn_src" "${pair_dir}/AD4Zn.dat" || return 1
    cp -f "$tz_src"    "${pair_dir}/${tz_file}" || return 1
    cd "$pair_dir" || return 1
    rm -f "$out_file"          # never leave a stale/partial result looking "done"

    # 4a: Grid parameter file (search box centred on this receptor's zinc)
    make_gpf "$tz_file" "$lig_abs" "$gpf_file" "$center" 30,30,30 || return 1

    # 4b: Generate affinity maps
    autogrid4 -p "$gpf_file" -l "$glg_file" || return 1

    # 4c: Clean up intermediates
    rm -f AD4Zn.dat "$tz_file"

    # 4d: Dock with Vina using AutoDock4 scoring
    vina --ligand "$lig_abs" --maps "${name}_tz" --scoring ad4 \
         --exhaustiveness 32 --cpu "$NCPU" \
         --out "$out_file" || { rm -f "$out_file"; return 1; }
}

# =========================================================
# STEP 3: PREPARE RECEPTORS  (skipped per-receptor if already done)
# =========================================================
# Supported formats: PDB or CIF
for receptor_file in receptors/*.{pdb,cif}; do
    base_name=$(basename "$receptor_file")
    name="${base_name%.*}"
    rec_dir="./receptors/${name}"
    rec_pdbqt="${rec_dir}/${name}.pdbqt"
    tz_pdbqt="${rec_dir}/${name}_tz.pdbqt"
    mkdir -p "$rec_dir"

    # 3a: Convert to PDBQT
    if up_to_date "$rec_pdbqt" "$receptor_file"; then
        echo "Skipping receptor conversion for $name (already done)"
    else
        echo "Preparing receptor: $receptor_file"
        rm -f "$rec_pdbqt"
        mk_prepare_receptor.py -i "$receptor_file" -o "${rec_dir}/${name}" -p \
            --allow_bad_res --default_altloc A \
            || { echo "mk_prepare_receptor.py failed for $name, skipping"; FAILED+=("receptor:$name"); continue; }
    fi

    # 3b: Add pseudoatoms for zinc coordination geometry
    if up_to_date "$tz_pdbqt" "$rec_pdbqt"; then
        echo "Skipping zinc pseudoatoms for $name (already done)"
    else
        rm -f "$tz_pdbqt"
        python "$VINA_SCRIPTS/zinc_pseudo.py" -r "$rec_pdbqt" -o "$tz_pdbqt" \
            || { echo "zinc_pseudo.py failed for $name, skipping"; FAILED+=("receptor:$name"); continue; }
    fi

    tz_src="$(realpath "$tz_pdbqt")"

    # 3c: Grid centre = first Zn atom. Read from the generated PDBQT (fixed
    #     PDB-style columns), which works for both PDB and CIF input.
    center=$(awk 'substr($0,1,6) ~ /^(ATOM|HETATM)/ && substr($0,13,4) ~ /ZN/ {
                      printf "%.3f,%.3f,%.3f\n", substr($0,31,8)+0, substr($0,39,8)+0, substr($0,47,8)+0
                      exit
                  }' "$rec_pdbqt")

    if [ -z "$center" ]; then
        echo "No Zn found in $receptor_file, skipping"
        FAILED+=("receptor:$name (no Zn)")
        continue
    fi
    echo "Grid centre for $name: $center"

    # -----------------------------------------------------
    # STEP 4: GENERATE GRID MAPS AND DOCK EACH LIGAND
    #         (a pair is skipped if its result already exists)
    # -----------------------------------------------------
    for lig in ligands/*/*.pdbqt; do
        lig_base=$(basename "$lig" .pdbqt)
        lig_abs="$(realpath "$lig")"
        pair_dir="${rec_dir}/${name}+${lig_base}"
        result="${pair_dir}/${name}_tz+${lig_base}.pdbqt"

        if up_to_date "$result" "$lig" "$tz_src"; then
            echo ">>> Skipping ${name} + ${lig_base} (result exists)"
            continue
        fi

        mkdir -p "$pair_dir"
        echo ">>> Processing ${name} + ${lig_base} in ${pair_dir}..."
        if ! ( dock_pair "$name" "$lig_abs" "$lig_base" "$tz_src" "$pair_dir" "$center" ); then
            echo "!!! Docking failed for ${name} + ${lig_base}"
            FAILED+=("pair:${name}+${lig_base}")
        fi
    done

    echo ">>> Receptor $name complete."
done

# =========================================================
# STEP 5: POST-PROCESSING
# =========================================================
if [ "${#FAILED[@]}" -gt 0 ]; then
    echo "The following steps failed:"
    printf '  %s\n' "${FAILED[@]}"
fi

cd "$JOB_DIR" || exit 1
if [ -z "$(find receptors -name '*_tz+*.pdbqt' -print -quit)" ]; then
    echo "No docking results found; skipping postprocess.py"
    exit 1
fi
python postprocess.py || { echo "postprocess.py failed"; exit 1; }

# Exit non-zero if anything failed, so SLURM records the problem
[ "${#FAILED[@]}" -eq 0 ]
