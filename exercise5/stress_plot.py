import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def plot_stenstor():
    filename = "stress_summary.txt"
    #elements are ratio , stress tensor11, stress tensor33
    data = np.loadtxt(filename, delimiter=",", skiprows=1)
    strain_rate = data[:,0]
    stress11 = data[:,1]
    stress33 = data[:,2]

    plt.figure()
    plt.scatter(strain_rate,stress11)
    plt.title("stress strain plot")
    plt.xlabel("strain rate")
    plt.ylabel("stress tensor 11")
    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.scatter(strain_rate,stress33)
    plt.title("stress strain plot")
    plt.xlabel("strain rate")
    plt.ylabel("stress tensor 33")
    plt.tight_layout()
    plt.show()


