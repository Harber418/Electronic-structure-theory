import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scienceplots

plt.style.use('science') # more scientific style for matplotlib
plt.rcParams['text.usetex'] = False # this avoids an annoying latex installation


def vinet_eos(v,v0,b0,b0prime,e0):
    nabla = (v/v0)**(1/3)
    return e0 + (4*b0*v0)/((b0prime-1)**2) + (2*b0*v0)/((b0prime-1)**2)*np.exp((3/2)*(b0prime-1)*(1-nabla))*(3*(b0prime-1)*(1-nabla)-2)


def fit(equ,energy,volume,type):

    v = volume
    e = energy
    B0 = 100*10**9 
    #100GPa
    #now convert to Ry/A^3 from joules per m^3
    #1 10^9 Pa = 10^9 J/m^3
    B0 = B0 * (10**(-30))/(1.60218*10**(-19)*13.60569312)
    #here we have the energy and volume of our initla sample 
    #in units of angstroms and Ry
    trial_v = 25


    #energy_trial = -19.19270380*13.605693 * 1.60218e-19
    #here is is in Ry
    trial_energy = -37
    
    #in units of angstrom
    p0 = [
    trial_v,         # v0
    B0,                     # b0 ~ Ry/A^3 
    4.0,                      # b0'
    trial_energy
]
    
    popt, _ = curve_fit(vinet_eos, v, e, p0=p0,maxfev=100000)
    

    v0, b0, b0p, e0 = popt

    print(f"paramertas for equation {equ}")
    print(f"E0  = {e0:.6f} Ry")
    print(f"V0  = {v0:.4f} Å^3")
    print(f"B0  = {b0:.5f} Ry/Å^3")
    print(f"B0' = {b0p:.4f}")
    vfit = np.linspace(v.min(), v.max(), 500)

    water_mass = (1.008*2 + 15.999)* 1.66054*10**(-27) #Kg
    titles="vinet"
    efit = vinet_eos(vfit, v0,b0,b0p,e0)
    vol_per_atom = v0
    print(f"the volume per atom is {vol_per_atom}")
    density = water_mass/(vol_per_atom*10**(-30))
    print(f"density is {density} in kg per m^3")
    B0 = b0 * (1.60218*10**(-19)*13.60569312)/(10**(-30))
    print(f"bulk is {B0*10**(-9)} in GPA")
    sound_velocity = np.sqrt(B0/density)
    print(f"the sound velocity is {sound_velocity} in m per s ")
    n = 1 / (vol_per_atom * 1e-30)
    kD = ((2*np.pi)**3*n*3/(4*np.pi))**(1/3)
    #hbar in units j s 
    #bolzmann constant in units m^2 kg s^-2 k^-1
    debye_temperature = (1.054571817*10**(-34)*sound_velocity*kD)/(1.380649 *10**(-23))
    print(f"the debye temperature is {debye_temperature} in kelvin ")

    colours = [False,"r","orange",False,False,False,False,False,"gold"]
    scatter = [False,"b","navy",False,False,False,False,False,"royalblue"]
    ices = [False,"Vinet : Ice I","Vinet : ice II",False,False,False,False,False,"Vinet : Ice VIII"]
    plt.figure()
    plt.scatter(volume, energy, s=25, label="DFT vdW",color = scatter[type])
    plt.plot(vfit, efit, 'r-', label=f"{ices[type]}",color=colours[type])
    plt.xlabel("Volume (Å³)")
    plt.ylabel("Energy (Ry)")
    plt.legend()
    plt.tight_layout()
    plt.show()

def joint_fit(equ,energy,volume,van_e,van_v,type):

    v = volume
    e = energy
    V = van_v
    E = van_e
    B0 = 18*10**9 
    #100GPa
    #now convert to Ry/A^3 from joules per m^3
    #1 10^9 Pa = 10^9 J/m^3
    B0 = B0 * (10**(-30))/(1.60218*10**(-19)*13.60569312)
    #here we have the energy and volume of our initla sample 
    #in units of angstroms and Ry
    trial_v = 25


    #energy_trial = -19.19270380*13.605693 * 1.60218e-19
    #here is is in Ry
    trial_energy = -37
    
    #in units of angstrom
    p0 = [
    trial_v,         # v0
    B0,                     # b0 ~ Ry/A^3 
    4.0,                      # b0'
    trial_energy
]
    
    popt, _ = curve_fit(vinet_eos, v, e, p0=p0,maxfev=100000)
    poptvan, _ = curve_fit(vinet_eos, V, E, p0=p0,maxfev=100000)

    v0, b0, b0p, e0 = popt
    v0van, b0van, b0pvan, e0van = poptvan
    print(f"paramertas for equation {equ}")
    print(f"E0  = {e0:.6f} Ry")
    print(f"V0  = {v0:.4f} Å^3")
    print(f"B0  = {b0:.5f} Ry/Å^3")
    print(f"B0' = {b0p:.4f}")
    print("=============================")

    print(f"paramertas for equation {equ}")
    print(f"E0  = {e0van:.6f} Ry")
    print(f"V0  = {v0van:.4f} Å^3")
    print(f"B0  = {b0van:.5f} Ry/Å^3")
    print(f"B0' = {b0pvan:.4f}")
    vfit = np.linspace(v.min(), v.max(), 500)

    si_mass = 28.085 * 1.66054*10**(-27) #Kg
    water_mass = (1.008*2 + 15.999)* 1.66054*10**(-27) #Kg
    titles="vinet"
    efit = vinet_eos(vfit, v0,b0,b0p,e0)
    vol_per_atom = v0
    print(f"the volume per atom is {vol_per_atom}")
    density = water_mass/(vol_per_atom*10**(-30))
    print(f"density is {density} in kg per m^3")
    B0 = b0 * (1.60218*10**(-19)*13.60569312)/(10**(-30))
    print(f"bulk is {B0*10**(-9)} in GPA")
    print("=============================")
    B0van = b0van * (1.60218*10**(-19)*13.60569312)/(10**(-30))
    print(f"bulk is {B0van*10**(-9)} in GPA for van")
    sound_velocity = np.sqrt(B0/density)
    print(f"the sound velocity is {sound_velocity} in m per s ")
    n = 1 / (vol_per_atom * 1e-30)
    kD = ((2*np.pi)**3*n*3/(4*np.pi))**(1/3)
    #hbar in units j s 
    #bolzmann constant in units m^2 kg s^-2 k^-1
    debye_temperature = (1.054571817*10**(-34)*sound_velocity*kD)/(1.380649 *10**(-23))
    print(f"the debye temperature is {debye_temperature} in kelvin ")

    colours = [False,"r","orange",False,False,False,False,"darkgoldenrod","gold"]
    scatter = [False,"b","navy",False,False,False,False,"midnightblue","royalblue"]
    ices = [False,"Vinet : Ice I","Vinet : ice II",False,False,False,False,False,"Vinet : Ice VIII"]
    Vfit = np.linspace(V.min(), V.max(), 500)

    Efit = vinet_eos(Vfit, v0van,b0van,b0pvan,e0van)

    plt.figure()
    plt.scatter(volume, energy, s=25, label="DFT PBS",color = scatter[type])
    plt.plot(vfit, efit, 'r-', label=f"{ices[type]} PBS",color=colours[type])

    plt.scatter(van_v, van_e, s=25, label="DFT VDW",color = scatter[type-1])
    plt.plot(Vfit, Efit, 'r-', label=f"{ices[type]} VDW",color=colours[type-1])
    plt.xlabel("Volume (Å³)", fontsize=16)
    plt.ylabel("Energy (Ry)", fontsize=16)
    plt.legend(fontsize=50)
    plt.tight_layout()
    plt.show()

def read_volume_energy(energy_file, volume_file):
    energies = []
    volumes = []
    # Read energies
    with open(energy_file, "r") as f:
        for line in f:
            #if "Final energy" in line:
            parts = line.split()
            energies.append(float(parts[4]))  # 5th element (index 4)
    # Read volumes
    with open(volume_file, "r") as f:
        for line in f:
            if "unit-cell volume" in line:
                parts = line.split()
                volume_au = float(parts[4])  # 5th element (index 4)
                volume_angstrom = volume_au * (0.529177 ** 3)
                volumes.append(volume_angstrom)
    energies = np.array(energies)
    volumes = np.array(volumes)

    return volumes, energies

def main():
    #vhange waht type of ice we measure 
    
    ice = 1




    v,e = read_volume_energy(f"energies_ice{ice}_van.txt",f"volume_ice{ice}_van.txt")
    #V,E = read_volume_energy(f"energies_ice{ice}_van.txt",f"volume_ice{ice}_van.txt")
    #v0_guess = v[np.argmin(e)]
    #mask = (v > 0.94*v0_guess) & (v < 1.06*v0_guess)

    #v = v[mask]
    #e = e[mask]
    #ice 1 
    #8 molecules 
    #ice 2 
    # 12 molecules 
    #ice 8 
    # 8 molecules 
    molecules =[0,8,12,0,0,0,0,0,8]
    v = v/molecules[ice]
    e = e/molecules[ice]
    #V = V/molecules[ice]
    #E= E/molecules[ice]
    #turn energy to j 
    #df["energy_J"] = df["energy_Ry"] * 13.605693 * 1.60218e-19


    #plt.figure()
    #plt.scatter(v, e, s=25)
    #plt.xlabel("Volume (Å³)")
    #plt.ylabel("Total energy (Ry)")
    #plt.tight_layout()
    #plt.show()
    #fit options 
    fit("vinet_eos",e, v,ice)
    #joint_fit("vinet_eos",e, v,E,V,ice)

if __name__ == "__main__":
    main()