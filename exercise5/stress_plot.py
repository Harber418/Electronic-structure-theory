import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def plot_stenstor():
    filename = "stress_summary.txt"
    #elements are ratio , stress tensor11, stress tensor33
    data = np.loadtxt(filename, delimiter=",", skiprows=1)
    strain_rate = data[:,0]*0.0001
    stress11 = data[:,1]
    stress33 = data[:,2]

    plt.figure()
    plt.scatter(strain_rate,stress11,c="b")
    plt.title("stress strain graph")
    plt.xlabel("strain rate")
    plt.ylabel("stress tensor 11")
    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.scatter(strain_rate,stress33,c="purple")
    plt.title("stress strain graph")
    plt.xlabel("strain rate")
    plt.ylabel("stress tensor 33")
    plt.tight_layout()
    plt.show()

def equation(x,m,c):
    return m*x +c


def fit():
    filename = "stress_summary.txt"
    #elements are ratio , stress tensor11, stress tensor33
    data = np.loadtxt(filename, delimiter=",", skiprows=1)
    strain_rate = data[:,0]*0.0001
    stress11 = data[:,1]
    stress33 = data[:,2]

    p0 = [1,0]
    popt, _ = curve_fit(equation, strain_rate, stress11, p0=p0,maxfev=100000)
    m,c = popt 

    y = equation(strain_rate,m,c)
    plt.figure()
    plt.scatter(strain_rate,stress11,c="b",s=2)
    plt.plot(strain_rate,y,c="r",)
    plt.title(f"stress strain graph. derivative = {m}")
    plt.xlabel("strain rate")
    plt.ylabel("stress tensor 11")
    plt.tight_layout()
    plt.show()


    popt, _ = curve_fit(equation, strain_rate, stress33, p0=p0,maxfev=100000)
    m,c = popt 

    y = equation(strain_rate,m,c)
    plt.figure()
    plt.scatter(strain_rate,stress33,c="purple",s=2)
    plt.plot(strain_rate,y,c="r")
    plt.title(f"stress strain graph. derivative = {m}")
    plt.xlabel("strain rate")
    plt.ylabel("stress tensor 33")
    plt.tight_layout()
    plt.show()



