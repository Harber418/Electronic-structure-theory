import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def plot_stenstor():
    filename = "stress.txt"
    #elements are ratio , stress tensor11, stress tensor33
    data = np.loadtxt(filename, delimiter=",", skiprows=1)
    #[;0] is 1 plus minis 1% of the tetragonal ratio 
    strain_rate = (data[:,0]*0.0001-1.0)/1.0
    #stress is in units Ry /borh ^3
    #should convert to gpa 
    stress11 = data[:,1]*14710.5
    stress33 = data[:,2]*14710.5
    #strain_rate = strain_rate[::-1]
    #stress11 = stress11[::-1]
    #stress33 = stress33[::-1]
    plt.figure()
    plt.scatter(strain_rate,stress11,c="b")
    plt.scatter(strain_rate,stress33,c="purple")
    plt.title("stress strain graph")
    plt.xlabel("strain rate")
    plt.ylabel("stress tensor")
    plt.tight_layout()
    plt.show()


def equation(x,m,c):
    return m*x +c


def fit():
    filename = "stress.txt"
    #elements are ratio , stress tensor11, stress tensor33
    data = np.loadtxt(filename, delimiter=",", skiprows=1)
    #ratio of lengths c/a = 1 +-1% where the length originally was 10.3335bhor 
    strain_rate = (data[:,0]*0.0001*10.3335 - 10.3335)/10.3335
    #units of Ry /bhr^3
    stress11 = data[:,1]*(1.60218*10**(-19)*13.60569312)
    stress33 = data[:,2]*(1.60218*10**(-19)*13.60569312)
    #14710.5
    p0 = [-1,0]
    popt, _ = curve_fit(equation, strain_rate, stress11, p0=p0,maxfev=100000)
    m,c = popt 

    y = equation(strain_rate,m,c)
    #-----------------------------
    popt2, _ = curve_fit(equation, strain_rate, stress33, p0=p0,maxfev=100000)
    m2,c2 = popt2 

    y2 = equation(strain_rate,m2,c2)
    plt.figure()
    plt.scatter(strain_rate,stress11,c="b",s=2,label="stress11")
    plt.plot(strain_rate,y,c="r",label=f"fit stress11. derivative = {m}")
    plt.scatter(strain_rate,stress33,c="purple",s=2,label="stress33")
    plt.plot(strain_rate,y2,c="orange",label=f"fit stress33. derivative = {m2}")
    plt.title(f"stress strain graph.")
    plt.xlabel("strain rate")
    plt.ylabel("stress tensor")
    plt.legend()
    plt.tight_layout()
    plt.show()



def main():
    plot_stenstor()
    fit()

if __name__ == "__main__":
    main()
