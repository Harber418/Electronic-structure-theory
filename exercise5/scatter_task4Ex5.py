import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
class fitting():
    def __init__(self):
        self.beta = []
        self.alpha = []
        
    def fit(self,equ,energy,volume,type ="si"):

        v = volume
        e = energy
        B0 = 100*10**9 
        #100GPa
        #now convert to Ry/A^3 from joules per m^3
        #1 10^9 Pa = 10^9 J/m^3
        B0 = B0 * (10**(-30))/(1.60218*10**(-19)*13.60569312)/2
        #here we have the energy and volume of our initla sample 
        #in units of angstroms and Ry
        trial_v = 20*2


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
        print(f"V0  = {v0/2:.4f} Å^3")
        print(f"B0  = {b0:.5f} Ry/Å^3")
        print(f"B0' = {b0p:.4f}")
        vfit = np.linspace(v.min(), v.max(), 500)

        si_mass = 28.085 * 1.66054*10**(-27) #Kg
        if type =="si":
            if equ == "murnagham_energy":
                efit = murnagham_energy(vfit, v0,b0,b0p,e0)
                titles="murnagham"
                vol_per_atom = v0/2
                print(f"the volume per atom is {vol_per_atom}")
                density = si_mass/(vol_per_atom*10**(-30))
                print(f"density is {density} in kg per m^3")
                B0 = b0 * (1.60218*10**(-19)*13.60569312)/(10**(-30))
                print(f"{B0*10**(-9)} in GPA")
                sound_velocity = np.sqrt(B0/density)
                print(f"the sound velocity is {sound_velocity} in m per s ")
                n = 1 / (vol_per_atom * 1e-30)
                kD = ((2*np.pi)**3*n*3/(4*np.pi))**(1/3)
                #hbar in units j s 
                #bolzmann constant in units m^2 kg s^-2 k^-1
                debye_temperature = (1.054571817*10**(-34)*sound_velocity*kD)/(1.380649 *10**(-23))
                print(f"the debye temperature is {debye_temperature} in kelvin ")
            elif equ == "BM_energy":
                efit = BM_energy(vfit, v0,b0,b0p,e0)
                titles="BM"
                vol_per_atom = v0/2
                print(f"the volume per atom is {vol_per_atom}")
                density = si_mass/(vol_per_atom*10**(-30))
                print(f"density is {density} in kg per m^3")
                B0 = b0 * (1.60218*10**(-19)*13.60569312)/(10**(-30))
                print(f"{B0*10**(-9)} in GPA")
                sound_velocity = np.sqrt(B0/density)
                print(f"the sound velocity is {sound_velocity} in m per s ")
                n = 1 / (vol_per_atom * 1e-30)
                kD = ((2*np.pi)**3*n*3/(4*np.pi))**(1/3)
                #hbar in units j s 
                #bolzmann constant in units m^2 kg s^-2 k^-1
                debye_temperature = (1.054571817*10**(-34)*sound_velocity*kD)/(1.380649 *10**(-23))
                print(f"the debye temperature is {debye_temperature} in kelvin ")
            elif equ == "vinet_eos":
                titles="vinet"
                efit = vinet_eos(vfit, v0,b0,b0p,e0)
                vol_per_atom = v0/2
                print(f"the volume per atom is {vol_per_atom}")
                density = si_mass/(vol_per_atom*10**(-30))
                print(f"density is {density} in kg per m^3")
                B0 = b0 * (1.60218*10**(-19)*13.60569312)/(10**(-30))
                print(f"{B0*10**(-9)} in GPA")
                sound_velocity = np.sqrt(B0/density)
                print(f"the sound velocity is {sound_velocity} in m per s ")
                n = 1 / (vol_per_atom * 1e-30)
                kD = ((2*np.pi)**3*n*3/(4*np.pi))**(1/3)
                #hbar in units j s 
                #bolzmann constant in units m^2 kg s^-2 k^-1
                debye_temperature = (1.054571817*10**(-34)*sound_velocity*kD)/(1.380649 *10**(-23))
                print(f"the debye temperature is {debye_temperature} in kelvin ")

        if type =="alpha":
            mass = 118.7* 1.66054*10**(-27)
            #alpha is diamond but has 2 atoms 
            titles="vinet"
            efit = vinet_eos(vfit, v0,b0,b0p,e0)
            vol_per_atom = v0/2
            print(f"the volume per atom is {vol_per_atom}")
            density = mass/(vol_per_atom*10**(-30))
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

        else:
            #beta     
            mass = 118.7* 1.66054*10**(-27)
            titles="vinet"
            efit = vinet_eos(vfit, v0,b0,b0p,e0)
            vol_per_atom = v0/2
            print(f"the volume per atom is {vol_per_atom}")
            density = mass/(vol_per_atom*10**(-30))
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

    

        plt.figure()
        plt.scatter(volume, energy, s=25, label="DFT")
        plt.plot(vfit, efit, 'r-', label=f"{titles} fit")
        plt.xlabel("Volume (Å³)")
        plt.ylabel("Energy (Ry)")
        plt.legend()
        plt.tight_layout()
        plt.show()

        return vfit ,efit,popt

    def read_volume_energy(self,energy_file, volume_file):
        energies = []
        volumes = []
        # Read energies
        with open(energy_file, "r") as f:
            for line in f:
                if "total energy" in line:
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

    def scatter(self,v,e):
        plt.figure()
        plt.scatter(v, e, s=25)
        plt.xlabel("Volume (Å³)")
        plt.ylabel("Total energy (Ry)")
        plt.tight_layout()
        plt.show()

    def all_together_now(self,v,e,Bv ,Be,vfita , efita,vfitb , efitb ):
        #we all live in a yellow submarine 
        E0 = np.min(efita)
        index = np.argmin(efita)
        Va = vfita[index]

        E0B = np.min(efitb)
        index = np.argmin(efitb)
        Vb = vfitb[index]
        #gradiant 
        m = np.round((E0-E0B)/(Va-Vb),3)

        plt.figure()
        plt.scatter(v, e, s=25,c="b")
        plt.scatter(Bv, Be, s=25,c="r")
        plt.plot(vfita,efita,label="alpha",c='purple')
        plt.plot(vfitb,efitb,label="beta",c="g")
        #plot the common tangent 
        plt.plot([Va,Vb],[E0,E0B],label =f"cotangent m = {m}",c='black')
        plt.legend()
        plt.tight_layout()
        plt.show()


    def A_derivative(self,v, v0,b0,b0p,e0):
        #calcualte the derivative at a given point 
        h = 1e-5
        right = v+h 
        left = v-h 
        return (Ealpha(right, v0,b0,b0p,e0)- Ealpha(left, v0,b0,b0p,e0))/(2*h)

    def B_derivative(self,v, v0,b0,b0p,e0):
        #calcualte the derivative at a given point 
        h = 1e-5
        right = v+h 
        left = v-h 
        return (Ebeta(right, v0,b0,b0p,e0)- Ebeta(left, v0,b0,b0p,e0))/(2*h)

    def transition(self, vars):
        #at the transition pressure the two phases have equal enthalpy
        # H = E +PV
        # so E + PV = E + PV 
        #p = Ea - Eb/va-vb
        #additionally
        #gradiants of the tangent at both points must be -P 
        Av,Bv,p = vars

        eq1 = self.A_derivative(Av, self.alpha[0],self.alpha[1],self.alpha[2],self.alpha[3]) + p
        eq2 = self.A_derivative(Bv, self.beta[0],self.beta[1],self.beta[2],self.beta[3]) + p

        eq3 = (Ebeta(Bv, self.beta[0],self.beta[1],self.beta[2],self.beta[3]) - Ealpha(Av, self.alpha[0],self.alpha[1],self.alpha[2],self.alpha[3]))/(Bv-Av) + p
        
        return [eq1, eq2, eq3]
    
    def solve(self,popta,poptb):
        initial_guess = [popta[0], poptb[0], 0.0]
        Va, Vb, P = fsolve(self.transition, initial_guess)
        return Va, Vb, P

def murnagham_energy(v,v0,b0,b0prime,e0):
        
        nabla = (v/v0)**(1/3)
        return  e0 + ((b0/b0prime)*v0*nabla**3)*(nabla**(-3*b0prime)/(b0prime -1)+1) -b0*v0/(b0prime-1)

def BM_energy(v,v0,b0,b0prime,e0):
    nabla = (v/v0)**(1/3)
    return e0 + (9/16)*b0*v0*(((nabla**(-2)-1)**3)*b0prime + (((nabla**(-2))-1)**2)*(6-4*(nabla**(-2))))


def vinet_eos(v,v0,b0,b0prime,e0):
    nabla = (v/v0)**(1/3)
    return e0 + (4*b0*v0)/((b0prime-1)**2) + (2*b0*v0)/((b0prime-1)**2)*np.exp((3/2)*(b0prime-1)*(1-nabla))*(3*(b0prime-1)*(1-nabla)-2)

def Ealpha(v,v0,b0,b0p,e0):
        #the alpha fit 
        return vinet_eos(v, v0,b0,b0p,e0)

def Ebeta(v,v0,b0,b0p,e0):
    #the beta fit 
    return vinet_eos(v, v0,b0,b0p,e0)

def main():

    run = fitting()

    #this is for alpha 
    v,e = run.read_volume_energy("energies_raw_alpha.txt","volume_data.alpha.txt")
    #but now for beta
    Bv ,Be = run.read_volume_energy("energies_raw_beta.txt","volume_data.beta.txt")
    #turn energy to j 
    #df["energy_J"] = df["energy_Ry"] * 13.605693 * 1.60218e-19


    run.scatter(v,e)
    vfita , efita ,popta =run.fit("vinet_eos",e, v,type ="alpha")
    v0, b0, b0p, e0 = popta
    run.alpha = [v0, b0, b0p, e0]
    run.scatter(Bv,Be)
    vfitb , efitb ,poptb= run.fit("vinet_eos",Be,Bv,type="beta")
    v0, b0, b0p, e0 = poptb
    run.beta = [v0, b0, b0p, e0]
    run.all_together_now(v,e,Bv ,Be,vfita , efita,vfitb , efitb)

    #I RAN NEW DATA FOR ALPHA BUT I NEED TO GET IT OFF THE SREVER
    Va, Vb, P = run.solve(popta,poptb)
    P_GPa = P * 14710.5
    print(P_GPa)
    deltaV = (Vb - Va)/Va * 100

    print(deltaV)




if __name__ == "__main__":
    main()
