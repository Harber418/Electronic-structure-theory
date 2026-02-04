import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
# -------------------------
# Read k-points
# -------------------------
k_data = []

with open("kpoints_totaltin.txt") as f:
    for line in f:
        parts = line.split()
        idx = int(parts[0].split("__")[1].split(".")[0])
        kpoints = int(parts[5])
        k_data.append([idx, kpoints])

df_k = pd.DataFrame(k_data, columns=["idx", "kpoints"])

# -------------------------
# Read energies
# -------------------------
e_data = []

with open("energy_total_tin.txt") as f:
    for line in f:
        parts = line.split()
        idx = int(parts[0].split("__")[1].split(".")[0])
        energy = float(parts[4])
        e_data.append([idx, energy])

df_e = pd.DataFrame(e_data, columns=["idx", "energy_Ry"])

# -------------------------
# Merge & sort
# -------------------------
df = (
    df_k
    .merge(df_e, on="idx")
    .sort_values("kpoints")
    .reset_index(drop=True)
)

# -------------------------
# Energy reference & ΔE
# -------------------------
E_ref = df["energy_Ry"].iloc[-1]
df["deltaE_meV"] = (df["energy_Ry"] - E_ref) * 13.605693 * 1000

# -------------------------
# Plot
# -------------------------
x = np.arange(1, 17)

plt.figure()
plt.plot(x, df["deltaE_meV"], marker="o")
plt.yscale("log")
plt.axhline(1, linestyle="--", color="r", label="1 meV")

plt.xlabel("Number of k-points")
plt.ylabel("Relative energy difference (meV)")
plt.title("Alpha-Tin k-point convergence")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(x, df["kpoints"], marker="o", c = "r", label="k points")
plt.plot(x,x**3, linestyle="--", label="cubic scaling", c="b")
plt.xlabel("Number of k-points")
plt.ylabel("K points")
plt.title("Alpha-Tin k-point convergence")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
