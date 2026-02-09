import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def murnagham_energy(v,v0,b0,b0prime,e0):
    
    nabla = (v/v0)**(1/3)
    return  e0 + ((b0/b0prime)*v0*nabla**3)*(nabla**(-3*b0prime)/(b0prime -1)+1) -b0*v0/(b0prime-1)

def BM_energy(v,v0,b0,b0prime,e0):
    nabla = (v/v0)**(1/3)
    return e0 + (9/16)*b0*v0*(((nabla**(-2)-1)**3)*b0prime + (((nabla**(-2))-1)**2)*(6-4*(nabla**-2)))


def vinet_eos(v,v0,b0,b0prime,e0):
    nabla = (v/v0)**(1/3)
    return e0 + (4*b0*v0)/((b0prime-1)**2) + (2*b0*v0)/((b0prime-1)**2)*np.exp((3/2)*(b0prime-1)*(1-nabla))*(3*(b0prime-1)*(1-nabla)-2)



def fit(equ,energy,volume):
    
    #mask = (volume > 900) & (volume < 1150)

    v = volume
    e = energy
    B0 = 100*10**9 
    #100GPa
    #now convert to Ry/A^3 from joules per m^3
    #1 10^9 Pa = 10^9 J/m^3
    B0 = B0 * (10**(-30))/(1.60218*10**(-19)*13.60569312)
    #here we have the energy and volume of our initla sample 
    #in units of angstroms and Ry
    trial_v = (543.09*10**(-2))**3

    #energy_trial = -19.19270380*13.605693 * 1.60218e-19
    #here is is in Ry
    trial_energy = -19.19270380
    #in units of angstrom
    p0 = [
    trial_v,         # v0
    B0,                     # b0 ~ Ry/A^3 
    4.0,                      # b0'
    trial_energy
]
    if equ == "murnagham_energy":
        popt, _ = curve_fit(murnagham_energy, v, e, p0=p0,maxfev=100000)
    elif equ == "BM_energy":
        popt, _ = curve_fit(BM_energy, v, e, p0=p0,maxfev=100000)
    elif equ == "vinet_eos":
        popt, _ = curve_fit(vinet_eos, v, e, p0=p0,maxfev=100000)
    

    v0, b0, b0p, e0 = popt

    print(f"paramertas for equation {equ}")
    print(f"E0  = {e0:.6f} Ry")
    print(f"V0  = {v0:.3f} Å^3")
    print(f"B0  = {b0:.4f} Ry/Å^3")
    print(f"B0' = {b0p:.3f}")
    vfit = np.linspace(v.min(), v.max(), 500)

    if equ == "murnagham_energy":
        efit = murnagham_energy(vfit, v0,b0,b0p,e0)
        titles="murnagham"
    elif equ == "BM_energy":
        efit = BM_energy(vfit, v0,b0,b0p,e0)
        titles="BM"
    elif equ == "vinet_eos":
        titles="vinet"
        efit = vinet_eos(vfit, v0,b0,b0p,e0)

    plt.figure()
    plt.scatter(volume, energy, s=25, label="DFT")
    plt.plot(vfit, efit, 'r-', label=f"{titles} fit")
    plt.xlabel("Volume (Å³)")
    plt.ylabel("Energy (Ry)")
    plt.legend()
    plt.tight_layout()
    plt.show()




def main():

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


    #turn to angstroms
    df["a_A"] = df["a"] * 0.529177
    df["volume_A"] = df["a_A"] **3
    #turn energy to j 
    #df["energy_J"] = df["energy_Ry"] * 13.605693 * 1.60218e-19


    plt.figure()
    plt.scatter(df["volume_A"], df["energy_Ry"], s=25)
    plt.xlabel("Volume (Å³)")
    plt.ylabel("Total energy (Ry)")
    plt.tight_layout()
    plt.show()
    #fit options 
    fit("murnagham_energy",df["energy_Ry"].values, df["volume_A"].values)
    fit("vinet_eos",df["energy_Ry"].values, df["volume_A"].values)
    fit("BM_energy",df["energy_Ry"].values, df["volume_A"].values)
if __name__ == "__main__":
    main()
