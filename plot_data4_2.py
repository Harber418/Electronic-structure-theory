import pandas as pd
import re
import matplotlib.pyplot as plt

alpha = True

if alpha:
    filename= "energies_raw_alpha_tin.txt"
    lines = 'alpha_tin_(\d+)\.out'
    title = "VAN"
    types = "alpha tin"
else:
    filename = "energies_raw_beta_tin.txt"
    lines = 'beta_tin__(\d+)\.out'
    title = "VAN"
    types = "beta tin"

data = []


with open(filename, "r") as f:
    for line in f:
        # Extract ECUT from filename
        ecut = int(re.search(lines, line).group(1))
        # Extract total energy
        energy = float(re.search(r'=\s*(-?\d+\.\d+)', line).group(1))
        data.append([ecut, energy])

data_beta = []
with open("energies_raw_beta_tin.txt", "r") as f:
    for line in f:
        # Extract ECUT from filename
        ecut = int(re.search('beta_tin__(\d+)\.out', line).group(1))
        # Extract total energy
        energy = float(re.search(r'=\s*(-?\d+\.\d+)', line).group(1))
        data_beta.append([ecut, energy])

# Create DataFrame
df = pd.DataFrame(data, columns=["ecut_Ry", "energy_Ry"])
df_beta = pd.DataFrame(data_beta, columns=["ecut_Ry", "energy_Ry"])

# Sort by cutoff
df = df.sort_values("ecut_Ry").reset_index(drop=True)
df_beta = df_beta.sort_values("ecut_Ry").reset_index(drop=True)

# Reference energy (highest cutoff)
#for alpha
E_ref = df["energy_Ry"].iloc[-1]
#for beta
E_ref_beta = df_beta["energy_Ry"].iloc[-1]

# Energy difference in meV 
df["deltaE_meV"] = (df["energy_Ry"] - E_ref) * 13.605693 * 1000
df_beta["deltaE_meV"] = (df_beta["energy_Ry"] - E_ref_beta) * 13.605693 * 1000

differnece = E_ref_beta - E_ref

# Plot
#make the y axis logarithmic
#add lines at 1 meV and 10 meV
plt.figure()
plt.plot(df["ecut_Ry"], df["deltaE_meV"], marker='o')
plt.yscale("log")
plt.axhline(5, linestyle='--', label="5 meV" , c='r')
#plt.axhline(10, linestyle='--', label="10 meV" , c='purple')
plt.xlabel("Plane-wave cutoff (Ry)")
plt.ylabel("Relative energy difference (meV)")
plt.title(f" {types} energy convergence, (potential : {title})")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

