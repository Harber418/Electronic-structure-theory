import numpy as np
import matplotlib.pyplot as plt
import re


def read_data(i):
    

    filename = f"quartz_dos_{i}.dat"   # change to beta, alpha, or quartz file

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
    
        plt.axvline(EF, linestyle='--')
    plt.xlabel("Energy (eV)")
    plt.ylabel("DOS")
    plt.title("DOS to Energy Beta-Tin")
    #plt.xlim(-10, 12)
    plt.legend()
    plt.show()


def get_fermi_from_scf(scf_file):
    with open(scf_file, "r") as f:
        for line in f:
            if "the Fermi energy is" in line:
                return float(line.split()[4])


def read_bands(filename, EF):
    k_points = []
    energies = []

    with open(filename, "r") as f:
        lines = f.readlines()

    # First line contains nbnd and nks
    header = lines[0]
    nbnd = int(header.split("nbnd=")[1].split(",")[0])
    
    i = 1
    while i < len(lines):
        # Read k-point
        k = list(map(float, lines[i].split()))
        k_points.append(k)
        i += 1

        # Read nbnd energies
        band_vals = []
        while len(band_vals) < nbnd:
            band_vals += list(map(float, lines[i].split()))
            i += 1

        energies.append(band_vals)

    k_points = np.array(k_points)
    energies = np.array(energies)

    # ---- Construct cumulative k-path distance ----
    k_dist = np.zeros(len(k_points))
    for j in range(1, len(k_points)):
        dk = np.linalg.norm(k_points[j] - k_points[j-1])
        k_dist[j] = k_dist[j-1] + dk

    # ---- Shift energies by Fermi energy ----
    #energies -= EF

    return k_dist, energies


def plot_bandsold():
    #EF = get_fermi_from_scf("quartz.scf.out")
    EF = 2.5756
    print("Fermi energy:", EF)

    k_dist, energies = read_bands("beta-bands.dat.gnu", EF)

    plt.figure()

    for band in range(energies.shape[1]):
        plt.plot(k_dist, energies[:, band])

    plt.axhline(EF, linestyle="--")
    plt.xlabel("k-path distance")
    plt.ylabel("Energy - EF (eV)")
    plt.title("Quartz Band Structure")
    plt.show()

def plot_bands(filename, EF):
    #EF = get_fermi_from_scf(scf_out)
    
    data = np.loadtxt(filename)

    k = data[:,0]
    energies = data[:,1:] 

    for i in range(energies.shape[1]):
        plt.plot(k, energies[:,i], color='black')

    plt.axhline(EF, linestyle='--')
    plt.xlabel("k-path")
    plt.ylabel("Energy (eV)")
    plt.title("Band Structure")
    plt.show()


def main():
    plots()
    #we need to put Fermi energy of each start 
    #alpha
    afe= 7.985
    filename = "beta-bands.dat.gru"
    #quartz
    Qfe= 5.230
    #beta 
    Bfe = 10.812
    plot_bands(filename,afe)


if __name__ == "__main__":
    main()