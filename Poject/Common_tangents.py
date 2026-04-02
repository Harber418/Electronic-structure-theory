import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import scienceplots

plt.style.use('science') # more scientific style for matplotlib
plt.rcParams['text.usetex'] = False # this avoids an annoying latex installation


class fitting():
    def __init__(self):
        self.ice1 = []
        self.ice2 = []

        self.ice8 = []
    def fit(self,equ,energy,volume,types):

        v = volume
        e = energy
        
        if types =="1":
            trial_v = 30
            B0 = 10*10**9 
            trial_energy = -34
            
        elif types =="2":
            trial_v = 30
            B0 = 10*10**9 
            trial_energy = -34

        elif types =="8":
            trial_v = 30
            B0 = 10*10**9 
            trial_energy = -34

        #now convert to Ry/A^3 from joules per m^3
        #1 10^9 Pa = 10^9 J/m^3
        B0 = B0 * (10**(-30))/(1.60218*10**(-19)*13.60569312)/2
        p0 = [
            trial_v,         # v0
            B0,                     # b0 ~ Ry/A^3 
            6.0,                      # b0'
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

        si_mass = 28.085 * 1.66054*10**(-27) #Kg

        equ == "vinet_eos"
        titles="vinet"
        efit = vinet_eos(vfit, v0,b0,b0p,e0)
        vol_per_atom = v0
        print(f"the volume  is {vol_per_atom} in angstrom cubed ")
        density = si_mass/(vol_per_atom*10**(-30))
        print(f"density is {density} in kg per m^3")
        B0 = b0 * ((1.60218*10**(-19))*13.60569312)/(10**(-30))
        print(f"{B0*10**(-9)} in GPA")
        sound_velocity = np.sqrt(B0/density)
        print(f"the sound velocity is {sound_velocity} in m per s ")
        n = 1 / (vol_per_atom * 1e-30)
        kD = ((2*np.pi)**3*n*3/(4*np.pi))**(1/3)
        #hbar in units j s 
        #bolzmann constant in units m^2 kg s^-2 k^-1
        debye_temperature = (1.054571817*10**(-34)*sound_velocity*kD)/(1.380649 *10**(-23))
        print(f"the debye temperature is {debye_temperature} in kelvin ")


        #plt.figure()
        #plt.scatter(volume, energy, s=25, label="DFT")
        #plt.plot(vfit, efit, 'r-', label=f"{titles} fit")
        #plt.xlabel("Volume (Å³)")
        #plt.ylabel("Energy (Ry)")
        #plt.legend()
        #plt.tight_layout()
        #plt.show()

        return vfit ,efit,popt

    def read_volume_energy(self,energy_file, volume_file):
        energies = []
        volumes = []
        # Read energies
        with open(energy_file, "r") as f:
            for line in f:
                #if "total energy" in line:
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

    def all_together_now(self, v, e, Bv, Be, v3, e3, vfita, efita, vfitb, efitb, vfit8, efit8,
                     Va, Vb, P, V2, V3, p2, types, Vone_to8, V8_to_one, p3):       #v1,e1,v2 ,e2,v3,e3,vfita , efita,vfitb , efitb,vfit8 , efit8,Va, Vb, P,V2, V3, p2,types
        #we all live in a yellow submarine 
        E0 = np.min(efita)
        index = np.argmin(efita)
        #Va = vfita[index]

        E0B = np.min(efitb)
        index = np.argmin(efitb)
        #Vb = vfitb[index]
        #gradiant 

        #line 1
        # Compute energies at transition volumes
        Ea = vinet_eos(Va, *self.ice1)
        Eb = vinet_eos(Vb, *self.ice2)
        E8 = vinet_eos(V3,*self.ice8)

        # Slope = -P
        slope = -P
        slope2 = -p2

        # Create volume range for plotting tangent
        Vline = np.linspace(min(Va, Vb)-5, max(Va, Vb)+5, 200)
        Vline2 = np.linspace(min(V2, V3)-5, max(V2, V3)+5, 200)
        Ea_1to8 = vinet_eos(Vone_to8, *self.ice1)
        E8_1to8 = vinet_eos(V8_to_one, *self.ice8)
        slope_1to8 = -p3
        Vline3 = np.linspace(min(Vone_to8, V8_to_one)-5, max(Vone_to8, V8_to_one)+5, 200)
        tangent_1to8_alpha = Ea_1to8 + slope_1to8 * (Vline3 - Vone_to8)
        tangent_1to8_beta  = E8_1to8 + slope_1to8 * (Vline3 - V8_to_one)        # Tangent lines
        tangent_alpha = Ea + slope * (Vline - Va)
        tangent_beta  = Eb + slope * (Vline - Vb)
        
        tangent_8 = E8 + slope2 *(Vline2-V3)
        
        colours = [False,"r","orange",False,False,False,False,False,"gold"]
        scatter = [False,"b","navy",False,False,False,False,False,"royalblue"]
        ices = [False,"Vinet : Ice I","Vinet : ice II",False,False,False,False,False,"Vinet : Ice VIII"]
        
        plt.figure()
        plt.scatter(v, e, s=15,color = scatter[types[0]],alpha=0.5)
        plt.scatter(Bv, Be, s=15,color = scatter[types[1]],alpha=0.5)
        plt.scatter(v3,e3,s=15,color = scatter[types[2]],alpha=0.5)
        plt.plot(vfita,efita,label=f"{ices[types[0]]}",color=colours[types[0]])
        plt.plot(vfitb,efitb,label=f"{ices[types[1]]}",color=colours[types[1]])
        plt.plot(vfit8,efit8,label=f"{ices[types[2]]}",color=colours[types[2]])
        plt.plot(Vline, tangent_alpha, 'k--', color = "k",label="Common Tangent 1↔2",alpha=0.8)
        plt.plot(Vline, tangent_beta, 'k--',color = "k",alpha= 0.9)
        plt.plot(Vline2, tangent_8, 'k-.', color = "gray",label="Common Tangent 2↔8",alpha=0.9)
        plt.plot(Vline3, tangent_1to8_alpha, 'g:', label="Common Tangent 1↔8", alpha=0.9)
        plt.plot(Vline3, tangent_1to8_beta, 'g:', alpha=0.9)
        plt.xlabel("Volume (Å³)")
        plt.tight_layout()
        plt.ylabel("Total energy (Ry)")
        
        #plot the common tangent 
        #plt.plot([Va,Vb],[E0,E0B],label =f"cotangent m = {m}",c='black')
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
        #p = -(Ea - Eb)/va-vb
        #additionally
        #gradiants of the tangent at both points must be -P 
        Av,Bv,p = vars

        eq1 = self.A_derivative(Av, self.ice1[0],self.ice1[1],self.ice1[2],self.ice1[3]) + p
        eq2 = self.B_derivative(Bv, self.ice2[0],self.ice2[1],self.ice2[2],self.ice2[3]) + p

        eq3 = (Ebeta(Bv, self.ice2[0],self.ice2[1],self.ice2[2],self.ice2[3]) - Ealpha(Av, self.ice1[0],self.ice1[1],self.ice1[2],self.ice1[3]))/(Bv-Av) + p
        
        return [eq1, eq2, eq3]
    
    def transition2to8(self, vars):
        #at the transition pressure the two phases have equal enthalpy
        # H = E +PV
        # so E + PV = E + PV 
        #p = -(Ea - Eb)/va-vb
        #additionally
        #gradiants of the tangent at both points must be -P 
        Av,Bv,p = vars

        eq1 = self.A_derivative(Av, self.ice2[0],self.ice2[1],self.ice2[2],self.ice2[3]) + p
        eq2 = self.B_derivative(Bv, self.ice8[0],self.ice8[1],self.ice8[2],self.ice8[3]) + p

        eq3 = (Ebeta(Bv, self.ice8[0],self.ice8[1],self.ice8[2],self.ice8[3]) - Ealpha(Av, self.ice2[0],self.ice2[1],self.ice2[2],self.ice2[3]))/(Bv-Av) + p
        
        return [eq1, eq2, eq3]
    
    def transition1to8(self, vars):
        #at the transition pressure the two phases have equal enthalpy
        # H = E +PV
        # so E + PV = E + PV 
        #p = -(Ea - Eb)/va-vb
        #additionally
        #gradiants of the tangent at both points must be -P 
        Av,Bv,p = vars

        eq1 = self.A_derivative(Av, self.ice1[0],self.ice1[1],self.ice1[2],self.ice1[3]) + p
        eq2 = self.B_derivative(Bv, self.ice8[0],self.ice8[1],self.ice8[2],self.ice8[3]) + p

        eq3 = (Ebeta(Bv, self.ice8[0],self.ice8[1],self.ice8[2],self.ice8[3]) - Ealpha(Av, self.ice1[0],self.ice1[1],self.ice1[2],self.ice1[3]))/(Bv-Av) + p
        
        return [eq1, eq2, eq3]
    
    def solve(self,popta,poptb,roll):
        initial_guess = [24,20,-0.01]
        if roll == "1":
            initial_guess = [popta[0] * 0.98, poptb[0] * 1.02, 0.0]
            Va, Vb, P = fsolve(self.transition, initial_guess)
        elif roll == "2":
            Va, Vb, P = fsolve(self.transition2to8, initial_guess)
        elif roll == "8":
            initial_guess = [24,20,-0.01]
            initial_guess = [popta[0] * 0.98, poptb[0] * 1.02, 0.0]          
            Va, Vb, P = fsolve(self.transition1to8, initial_guess)
        return Va, Vb, P


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
    #ice 1 
    #8 molecules 
    #ice 2 
    # 12 molecules 
    #ice 8 
    # 8 molecules 
    run = fitting()
    molecules =[0,8,12,0,0,0,0,0,8]
    ice1 = 1
    v1,e1 = run.read_volume_energy(f"energies_ice{ice1}.txt",f"volume_ice{ice1}.txt")
    v1 = v1/molecules[ice1]
    e1 = e1/molecules[ice1]
    ice2 = 2
    v2,e2 = run.read_volume_energy(f"energies_ice{ice2}.txt",f"volume_ice{ice2}.txt")
    v2 = v2/molecules[ice2]
    e2 = e2/molecules[ice2]

    ice3 = 8
    v3,e3 = run.read_volume_energy(f"energies_ice{ice3}.txt",f"volume_ice{ice3}.txt")
    v3 = v3/molecules[ice3]
    e3 = e3/molecules[ice3]
    #run.scatter(v,e)
    print("we now run ice 1 =========================")
    vfita , efita ,poptice1 =run.fit("vinet_eos",e1, v1,"1")
    v0, b0, b0p, e0 = poptice1
    run.ice1 = [v0, b0, b0p, e0]
    #run.scatter(Bv,Be)
    print("we now run ice 2 =========================")

    vfitb , efitb ,poptice2= run.fit("vinet_eos",e2,v2,"2")
    v0, b0, b0p, e0 = poptice2
    run.ice2 = [v0, b0, b0p, e0]
    #run for ice 8 
    print("we now run ice 8 =========================")

    vfit8 , efit8 ,poptice8= run.fit("vinet_eos",e3,v3,"8")
    v0, b0, b0p, e0 = poptice8
    run.ice8 = [v0, b0, b0p, e0]

    types =[ice1,ice2,ice3]
    #-------------------------
    #transition for ICE1 to ICE 2 
    #------------------------------
    Va, Vb, P = run.solve(poptice1,poptice2,"1")

    #========================
    #now from ice 2 to ice 8 
    #========================
    V2, V3, p2 = run.solve(poptice2,poptice8,"2")
    #========================
    #now from ice 1 to ice 8 
    #========================
    Vone_to8, V8_to_one, p3 = run.solve(poptice1,poptice8,"8")
    #================
    #and plot together 
    #====================



    run.all_together_now(
        v1, e1, v2, e2, v3, e3,
        vfita, efita, vfitb, efitb, vfit8, efit8,
        Va, Vb, P, V2, V3, p2, types,
        Vone_to8, V8_to_one, p3
    )
    print("solving ")
    #P_GPa = P * 14710.5
    P_GPa = P * (13.605) *(1.602e-19)/1e-30/1e9
    print("ice 1 to 2 ")
    print(f"transition pressure is {P_GPa}")
    print("==================================")
    P_GPa2 = p2 * (13.605) *(1.602e-19)/1e-30/1e9
    print("ice 2 to 8 ")
    print(f"transition pressure is {P_GPa2}")
    print("==================================")
    P_GPa2 = p3 * (13.605) *(1.602e-19)/1e-30/1e9
    print("ice 1 to 8 ")
    print(f"transition pressure is {P_GPa2}")


    print("va is ")
    print(Va)
    print(Vb)
    deltaV = (Va - Vb)/Va * 100

    print(f" delta v is {deltaV}")
    

if __name__ == "__main__":
    main()
