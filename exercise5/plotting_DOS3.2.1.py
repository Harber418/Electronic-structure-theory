import numpy as np
import matplotlib.pyplot as plt
import re


def read_data(i):
    

    filename = f"beta_sn_dos_{i}.dat"   # change to beta, alpha, or quartz file

    # ---- Read Fermi energy from header ----
    with open(filename, "r") as f:
        first_line = f.readline()

    # Extract number after "EFermi ="
    EF = float(re.search(r'EFermi\s*=\s*([0-9.]+)', first_line).group(1))

    print(f"Fermi energy extracted: {EF} eV")

    # ---- Load numerical data ----
    data = np.loadtxt(filename, comments="#")

    E = data[:,0]
    DOS = data[:,1]

    # ---- Shift energy axis ----
    E_shift = E - EF

    return E_shift, DOS,EF


def plots():
    E_data=[]
    DOS_data=[]
    # ---- Plot ----
    plt.figure()
    #4,6,8,10,12,14,16,18,20
    for i in range(4,22,2):
        E , DOS,EF = read_data(i)
        E_data.append(E)
        DOS_data.append(DOS)
        plt.plot(E, DOS,label=f"k point grid size: {i}")
    #plt.plot(E_data[1], DOS_data[1],label="k point grid size: 6")
    #plt.plot(E_data[2], DOS_data[2],label="k point grid size: 8")
    #plt.plot(E_data[3], DOS_data[3],label="k point grid size: 10")
  #  plt.plot(E_data[4], DOS_data[4],label="k point grid size: 12")
   # plt.plot(E_data[6], DOS_data[6],label="k point grid size: 16")
    #plt.plot(E_data[7], DOS_data[7],label="k point grid size: 18")
   # plt.plot(E_data[5],DOS_data[5],label ="k point grid size: 14")
   # plt.plot(E_data[8],DOS_data[8],label ="k point grid size: 20")
        plt.axvline(EF, linestyle='--')
    plt.xlabel("Energy (eV)")
    plt.ylabel("DOS")
    plt.title("DOS to Energy Beta-Tin")
    #plt.xlim(-10, 12)
    plt.legend()
    plt.show()

def main():
    plots()


if __name__ == "__main__":
    main()