import os, glob, csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ===============================
#  Utilities
# ===============================y

def get_top5_energies(pdbqt_path):
    """Extract top 5 docking energies from a PDBQT file."""
    energies = []
    try:
        with open(pdbqt_path) as f:
            for line in f:
                if line.strip().startswith("REMARK VINA RESULT:"):
                    try:
                        energies.append(float(line.split()[3]))
                    except Exception:
                        pass
        energies.sort()
        return energies[:5] if energies else None
    except Exception:
        return None


def get_atoms(pdb_path):
    """Extract Zn and ligand atom coordinates (only O/N/S for ligand)."""
    zn_atoms, lig_atoms = [], []
    protein_residues = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
        "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
        "THR","TRP","TYR","VAL"
    }
    allowed_elems = {"O", "N", "S"}

    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                elem = atom_name[0].upper()

                if res_name.upper() == "ZN" or "ZN" in atom_name.upper():
                    zn_atoms.append(np.array([x, y, z]))
                elif res_name.upper() not in protein_residues and elem in allowed_elems:
                    lig_atoms.append(np.array([x, y, z]))

    return zn_atoms, lig_atoms


def compute_min_distance(pdb_path):
    """Compute minimum Zn–ligand (O/N/S only) distance."""
    zn_coords, lig_coords = get_atoms(pdb_path)
    if not zn_coords or not lig_coords:
        return None
    dists = [np.linalg.norm(z - l) for z in zn_coords for l in lig_coords]
    return float(min(dists)) if dists else None


def merge_pdbs(receptor_file, ligand_file, output_file):
    """Merge receptor and ligand coordinates into one PDB."""
    with open(output_file, "w") as out:
        for src in [receptor_file, ligand_file]:
            with open(src) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        out.write(line[:80] + "\n")
        out.write("END\n")


# ===============================
#  Main analysis
# ===============================

os.makedirs("analysis_outputs/merged_structures", exist_ok=True)
results = []

for rec_dir in sorted(glob.glob("receptors/*/")):
    receptor_pdbqt_list = glob.glob(os.path.join(rec_dir, "*_tz.pdbqt"))
    if not receptor_pdbqt_list:
        print(f"[!] No receptor file found in {rec_dir}")
        continue

    receptor_pdbqt = receptor_pdbqt_list[0]
    receptor_name = os.path.basename(os.path.normpath(rec_dir))

    for lig_dir in sorted(glob.glob(os.path.join(rec_dir, "*+*/"))):
        ligand_pdbqt_list = glob.glob(os.path.join(lig_dir, "*.pdbqt"))
        if not ligand_pdbqt_list:
            continue

        ligand_pdbqt = ligand_pdbqt_list[0]
        ligand_name = os.path.basename(os.path.normpath(lig_dir)).split("+")[-1]

        energies = get_top5_energies(ligand_pdbqt)
        if not energies:
            print(f"[!] No energies found for {ligand_name}")
            continue

        energy_mean = np.mean(energies)
        energy_var = np.std(energies)

        merged_name = f"{receptor_name}_{ligand_name}_merged.pdb"
        merged_path = os.path.join("analysis_outputs/merged_structures", merged_name)
        merge_pdbs(receptor_pdbqt, ligand_pdbqt, merged_path)

        dist = compute_min_distance(merged_path)
        if dist is None:
            print(f"[!] Could not compute Zn-ligand distance for {receptor_name}+{ligand_name}")

        results.append({
            "pair": f"{receptor_name}+{ligand_name}",
            "receptor": receptor_name,
            "ligand": ligand_name,
            "top5_mean_energy": energy_mean,
            "top5_var_energy": energy_var,
            "min_Zn_ligand_distance_A": dist
        })


# ===============================
#  Aggregate across receptors
# ===============================

agg_by_ligand = {}
for r in results:
    lig = r["ligand"]
    agg_by_ligand.setdefault(lig, []).append(r["top5_mean_energy"])

for lig, vals in agg_by_ligand.items():
    mean_between = np.mean(vals)
    var_between = np.std(vals)
    for r in results:
        if r["ligand"] == lig:
            r["mean_between_receptors"] = mean_between
            r["var_between_receptors"] = var_between


# ===============================
#  Output results
# ===============================

csv_path = "analysis_outputs/docking_summary.csv"
os.makedirs("analysis_outputs", exist_ok=True)

with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "pair","receptor","ligand",
        "top5_mean_energy","top5_var_energy",
        "mean_between_receptors","var_between_receptors",
        "min_Zn_ligand_distance_A"
    ])
    writer.writeheader()
    writer.writerows(results)

print(f"✅ Results written to '{csv_path}'")


# ===============================================
# Overlapping receptor comparison: tallest-at-back plotting
# ===============================================
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os

# === Load data ===
csv_path = "analysis_outputs/docking_summary.csv"
df = pd.read_csv(csv_path)
df["top5_mean_energy"] = pd.to_numeric(df["top5_mean_energy"], errors="coerce")

# --- Identify substrate/hydrolysed pairs ---
substrates = [lig for lig in df["ligand"].unique() if not lig.startswith("hydrolysed_")]
substrates_mean = df[df["ligand"].isin(substrates)].groupby("ligand")["top5_mean_energy"].mean()
substrates_sorted = substrates_mean.sort_values(ascending=False).index.tolist()

# Build final ligand order: substrate first, hydrolysed second
ligands_sorted = []
for s in substrates_sorted:
    ligands_sorted.append(s)
    hydro = f"hydrolysed_{s}"
    if hydro in df["ligand"].unique():
        ligands_sorted.append(hydro)

# --- Receptors ---
receptors = sorted(df["receptor"].unique())

# --- Plot grouped bar chart ---
plt.figure(figsize=(12, 6))
bar_width = 0.12
x = np.arange(len(ligands_sorted))

for i, receptor in enumerate(receptors):
    subset = df[df["receptor"] == receptor].set_index("ligand").reindex(ligands_sorted)
    plt.bar(
        x + i * bar_width,
        subset["top5_mean_energy"],
        width=bar_width,
        label=receptor,
        alpha=0.8,
        edgecolor="black",
    )

plt.xticks(x + bar_width * (len(receptors)-1)/2, ligands_sorted, rotation=45, ha="right")
plt.ylabel("Mean Docking Energy (kcal/mol)")
plt.title("Binding Energies Across Receptors (Grouped by Ligand)")
plt.legend(title="Receptor", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()

os.makedirs("analysis_outputs", exist_ok=True)
out_path = "analysis_outputs/grouped_energy_per_ligand_sorted.png"
plt.savefig(out_path, dpi=300)
plt.close()
print(f"✅ Grouped bar chart saved as '{out_path}'")

# Pivot table for heatmap
pivot = df.pivot_table(index="receptor", columns="ligand", values="top5_mean_energy", aggfunc="mean")
pivot = pivot.reindex(columns=ligands_sorted)
pivot = pivot.replace(np.nan, -10)

pivot = pivot[~pivot.index.str.endswith("_1metal")]
pivot = pivot.loc[:, ~pivot.columns.str.startswith("hydrolysed_")]


plt.figure(figsize=(14, 6))
sns.heatmap(
    pivot,
    cmap="coolwarm_r",
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    cbar_kws={"label": "Mean Binding Energy (kcal/mol)"}
)
plt.title("Docking Energies: Receptor vs Ligand (Sorted by Substrate Energy)")
plt.xlabel("Ligand")
plt.ylabel("Receptor")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

out_path2 = "analysis_outputs/heatmap_energy_per_receptor_ligand_sorted.png"
plt.savefig(out_path2, dpi=300)
plt.close()
print(f"✅ Heatmap saved as '{out_path2}'")





plt.figure()

out_path3 = "analysis_outputs/clustermap_energy_per_receptor_ligand_sorted.png"



g = sns.clustermap(
    pivot,
    cmap="coolwarm_r",
    annot=True,
    method="average",
    #fmt=".2f",
    #linewidths=0.5,
    cbar_kws={"label": "Mean Binding Energy (kcal/mol)"}
)

g.fig.suptitle("Docking Energies: Receptor vs Ligand (Sorted by Substrate Energy)")

g.savefig(out_path3, dpi=300)
plt.close(g.fig)
