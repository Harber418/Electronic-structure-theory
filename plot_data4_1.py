import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np 
mt = True

if mt:
    filename= "energies_raw.txt"
    lines = 'quartz_mt_(\d+)\.out'
    title = "martins-troullier"
else:
    filename = "energies_raw_us.txt"
    lines = 'quartz_us_(\d+)\.out'
    title = "ultra soft"



data = []


with open(filename, "r") as f:
    for line in f:
        # Extract ECUT from filename
        ecut = int(re.search(lines, line).group(1))
        # Extract total energy
        energy = float(re.search(r'=\s*(-?\d+\.\d+)', line).group(1))
        data.append([ecut, energy])

# Create DataFrame
df = pd.DataFrame(data, columns=["ecut_Ry", "energy_Ry"])

# Sort by cutoff
df = df.sort_values("ecut_Ry").reset_index(drop=True)

# Reference energy (highest cutoff)
E_ref = df["energy_Ry"].iloc[-1]

# Energy difference in meV per SiO2
df["deltaE_meV"] = (df["energy_Ry"] - E_ref) * 13.605693 * 1000

# Plot
#make the y axis logarithmic
#add lines at 1 meV and 10 meV
plt.figure()
plt.plot(df["ecut_Ry"], df["deltaE_meV"], marker='o')
plt.yscale("log")
plt.axhline(1, linestyle='--', label="1 meV" , c='r')
plt.axhline(10, linestyle='--', label="10 meV" , c='purple')
plt.xlabel("Plane-wave cutoff (Ry)")
plt.ylabel("Relative energy difference (meV)")
plt.title(f"Quartz (SiO2) energy convergence potential : {title}")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
