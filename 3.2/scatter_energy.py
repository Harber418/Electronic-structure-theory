import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

filename = "energies_raw_us.txt"
pattern = r'si_volume_(\d+)\.out'

data = []

with open(filename) as f:
    for line in f:
        m1 = re.search(pattern, line)
        m2 = re.search(r'=\s*(-?\d+\.\d+)', line)
        if m1 and m2:
            radius = int(m1.group(1)) * 0.01
            energy = float(m2.group(1))
            data.append([radius, energy])

df = pd.DataFrame(data, columns=["a", "energy_Ry"])
df = df.sort_values("a").reset_index(drop=True)

df["volume"] = df["a"]**3
#turn to angstroms
df["a_A"] = df["a"] * 0.529177
df["volume_A"] = df["a_A"] **3
#turn energy to j 
df["energy_J"] = df["energy_Ry"] * 13.605693 * 1.60218e-19


plt.figure()
plt.scatter(df["volume_A"], df["energy_J"], s=25)
plt.xlabel("Volume (Å³)")
plt.ylabel("Total energy (J)")
plt.tight_layout()
plt.show()

def murnagham_energy(v,v0,b0,b0prime,e0):
    
    nabla = (v/v0)**(1/3)
    return  e0 + ((b0/b0prime)*v0*nabla**3)*(nabla**(-3*b0prime)/(b0prime -1)+1) -b0*v0/(b0prime-1)

def BM_energy(v,v0,b0,b0prime,e0):
    nabla = (v/v0)**(1/3)
    return e0 + (9/16)*b0*v0*(((nabla**(-2)-1)**3)*b0prime + (((nabla**(-2))-1)**2)*(6-4*nabla**-2))



def vinet_eos(v,v0,b0,b0prime,e0):
    nabla = (v/v0)**(1/3)
    return e0 + (4*b0*v0)/((b0prime-1)**2) + (2*b0*v0)/((b0prime-1)**2)*np.exp((3/2)*(b0prime-1)*(1-nabla))*(3*(b0prime-1)*(1-nabla)-2)



def fit(equ):
    mask = (df["volume"] > 900) & (df["volume"] < 1150)

    v = df.loc[mask, "volume_A"].values
    e = df.loc[mask, "energy_J"].values
    B0 = 100*10**9 
    #100GPa
    trial_v = (543.09*10**(-2))**3
    energy_trial = -19.19270380*13.605693 * 1.60218e-19
    #in units of angstrom
    p0 = [
    energy_trial,                 # e0
    trial_v,         # v0
    B0,                     # b0 ~ Ry/Bohr^3 (reasonable)
    4.0                      # b0'
]
    popt, pcov = curve_fit(equ, v, e, p0=p0)

    e0, v0, b0, b0p = popt

    print(f"E0  = {e0:.6f} Ry")
    print(f"V0  = {v0:.3f} Bohr^3")
    print(f"B0  = {b0:.4f} Ry/Bohr^3")
    print(f"B0' = {b0p:.3f}")
    vfit = np.linspace(v.min(), v.max(), 500)
    efit = equ(vfit, *popt)

    plt.figure()
    plt.scatter(df["volume"], df["energy_Ry"], s=25, label="DFT")
    plt.plot(vfit, efit, 'r-', label="BM fit")
    plt.xlabel("Volume (Bohr³)")
    plt.ylabel("Energy (Ry)")
    plt.legend()
    plt.tight_layout()
    plt.show()




def main():
    fit(BM_energy)

if __name__ == "__main__":
    main()
