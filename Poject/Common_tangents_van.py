import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, fsolve

def vinet_eos(v, v0, b0, b0prime, e0):
    nabla = (v / v0) ** (1 / 3)
    return e0 + (4 * b0 * v0) / ((b0prime - 1) ** 2) + (2 * b0 * v0) / ((b0prime - 1) ** 2) * np.exp((3 / 2) * (b0prime - 1) * (1 - nabla)) * (3 * (b0prime - 1) * (1 - nabla) - 2)

def read_volume_energy(energy_file, volume_file, molecules):
    energies = []
    volumes = []
    with open(energy_file, "r") as f:
        for line in f:
            if "total energy" in line:
                parts = line.split()
                energies.append(float(parts[4]))
    with open(volume_file, "r") as f:
        for line in f:
            if "unit-cell volume" in line:
                parts = line.split()
                volume_au = float(parts[4])
                volume_angstrom = volume_au * (0.529177 ** 3)
                volumes.append(volume_angstrom)
    energies = np.array(energies) / molecules
    volumes = np.array(volumes) / molecules
    return volumes, energies

def fit_eos(volumes, energies):
    B0 = 100 * 10 ** 9
    B0 = B0 * (10 ** (-30)) / (1.60218 * 10 ** (-19) * 13.60569312) 
    trial_v = np.mean(volumes)
    trial_energy = np.min(energies)
    p0 = [trial_v, B0, 4.0, trial_energy]
    popt, _ = curve_fit(vinet_eos, volumes, energies, p0=p0, maxfev=100000)
    v0, b0, b0p, e0 = popt
    print("==============")
    print(f"paramertas for equation ")
    print(f"E0  = {e0:.6f} Ry")
    print(f"V0  = {v0:.4f} Å^3")
    print(f"B0  = {b0:.5f} Ry/Å^3")
    print(f"B0' = {b0p:.4f}")
    B0 = b0 * (1.60218*10**(-19)*13.60569312)/(10**(-30))
    print(f"{B0*10**(-9)} in GPA")
    vfit = np.linspace(volumes.min(), volumes.max(), 500)
    efit = vinet_eos(vfit, *popt)
    return vfit, efit, popt

def tangent_line(v, v0, b0, b0p, e0, Vt, slope):
    Et = vinet_eos(Vt, v0, b0, b0p, e0)
    return Et + slope * (v - Vt)

def eos_derivative(v, v0, b0, b0p, e0):
    h = 1e-5
    return (vinet_eos(v + h, v0, b0, b0p, e0) - vinet_eos(v - h, v0, b0, b0p, e0)) / (2 * h)

def transition_eq(vars, poptA, poptB):
    Av, Bv, p = vars
    eq1 = eos_derivative(Av, *poptA) + p
    eq2 = eos_derivative(Bv, *poptB) + p
    eq3 = (vinet_eos(Bv, *poptB) - vinet_eos(Av, *poptA)) / (Bv - Av) + p
    return [eq1, eq2, eq3]

def solve_transition(poptA, poptB):
    initial_guess = [poptA[0], poptB[0], 0.0]
    Va, Vb, P = fsolve(transition_eq, initial_guess, args=(poptA, poptB))
    return Va, Vb, P

def plot_all(phases, fits, transitions, labels, colors, linestyles):
    plt.figure(figsize=(10, 7))
    # Scatter and fits
    for i, (v, e) in enumerate(phases):
        plt.scatter(v, e, color=colors[i], label=labels[i], alpha=0.7, marker='o')
        plt.plot(fits[i][0], fits[i][1], color=colors[i], linestyle=linestyles[i], label=labels[i] + " fit")
    # Tangents
    for idx, (V1, V2, P, poptA, poptB, color, label) in enumerate(transitions):
        vline = np.linspace(min(V1, V2) - 2, max(V1, V2) + 2, 100)
        plt.plot(vline, tangent_line(vline, *poptA, V1, -P), color=color, linestyle='--', label=label + " tangent A")
        plt.plot(vline, tangent_line(vline, *poptB, V2, -P), color=color, linestyle=':', label=label + " tangent B")
    plt.xlabel("Volume (Å³)")
    plt.ylabel("Total energy (Ry)")
    plt.legend()
    plt.tight_layout()
    plt.title("EoS fits and common tangents for normal and van data")
    plt.show()

def main():
    # File and molecule info
    phases = [
        ("energies_ice1.txt", "volume_ice1.txt", 8, "Ice I normal"),
        ("energies_ice2.txt", "volume_ice2.txt", 12, "Ice II normal"),
        ("energies_ice8.txt", "volume_ice8.txt", 8, "Ice VIII normal"),
        ("energies_ice1_van.txt", "volume_ice1_van.txt", 8, "Ice I van"),
        ("energies_ice2_van.txt", "volume_ice2_van.txt", 12, "Ice II van"),
        ("energies_ice8_van.txt", "volume_ice8_van.txt", 8, "Ice VIII van"),
    ]
    colors = ["red", "orange", "gold", "blue", "green", "purple"]
    linestyles = ["-", "-", "-", "--", "--", "--"]
    # Read and fit all phases
    phase_data = []
    fits = []
    popts = []
    labels = []
    for i, (e_file, v_file, n_mol, label) in enumerate(phases):
        print(f"ice {i} 3 normal 3 van ")
        v, e = read_volume_energy(e_file, v_file, n_mol)
        phase_data.append((v, e))
        vfit, efit, popt = fit_eos(v, e)
        fits.append((vfit, efit))
        popts.append(popt)
        labels.append(label)
    # Compute transitions (normal: 1-2, 2-8; van: 1-2, 2-8)
    transitions = []
    # Normal
    Va, Vb, P = solve_transition(popts[0], popts[1])
    transitions.append((Va, Vb, P, popts[0], popts[1], "black", "Normal 1→2"))
    Va, Vb, P = solve_transition(popts[1], popts[2])
    transitions.append((Va, Vb, P, popts[1], popts[2], "grey", "Normal 2→8"))
    # Van
    Va, Vb, P = solve_transition(popts[3], popts[4])
    transitions.append((Va, Vb, P, popts[3], popts[4], "navy", "Van 1→2"))
    Va, Vb, P = solve_transition(popts[4], popts[5])
    transitions.append((Va, Vb, P, popts[4], popts[5], "teal", "Van 2→8"))
    # Plot all
    plot_all(phase_data, fits, transitions, labels, colors, linestyles)

if __name__ == "__main__":
    main()