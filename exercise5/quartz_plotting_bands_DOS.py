import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib import gridspec

plt.rcParams["figure.dpi"] = 150
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["figure.figsize"] = (8, 6)

# ==========================
# ---- BAND STRUCTURE ----
# ==========================

band_data = np.loadtxt('quartz-bands.dat.gnu')

k = np.unique(band_data[:,0])
nk = len(k)

total_points = len(band_data)
nbands = total_points // nk

bands = band_data[:,1][:nbands*nk].reshape((nbands, nk))

# ==========================
# -------- DOS ------------
# ==========================

filename = "quartz_dos_12.dat"

with open(filename, "r") as f:
    first_line = f.readline()

EF = float(re.search(r'EFermi\s*=\s*([0-9.]+)', first_line).group(1))

dos_data = np.loadtxt(filename, comments="#")
E = dos_data[:,0]
DOS = dos_data[:,1]

# ==========================
# -------- PLOTTING -------
# ==========================

fig = plt.figure()
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05)

ax_band = plt.subplot(gs[0])
ax_dos  = plt.subplot(gs[1], sharey=ax_band)

# ---- Plot bands ----
for band in range(nbands):
    ax_band.plot(k, bands[band,:], linewidth=1, alpha=0.5, color='k')

# Fermi energy line
ax_band.axhline(EF, linestyle=(0, (5, 5)), linewidth=0.75, color='k', alpha=0.7)
ax_band.text(k[-1]*0.7, EF+0.2, r"$E_F$", fontsize=9)

# High symmetry lines
sym_points = [0, 0.5774, 0.9107, 1.5773, 2.0319, 2.6092, 2.9426, 3.7494]
for x in sym_points:
    ax_band.axvline(x, linewidth=0.75, color='k', alpha=0.5)

ax_band.set_xticks(sym_points)
ax_band.set_xticklabels(['$\Gamma$', 'M', 'K', '$\Gamma$', 'A', 'L', 'H', '$\Gamma$'])
ax_band.set_ylabel("Energy (eV)")
ax_band.set_xlim(min(k), max(k))


# ---- Plot DOS ----
ax_dos.plot(DOS, E, color='k', linewidth=1)
ax_dos.axhline(EF, linestyle=(0, (5, 5)), linewidth=0.75, color='k', alpha=0.7)

ax_dos.set_xlabel("DOS")
ax_dos.set_xlim(left=0)

# Remove duplicate y-axis labels
plt.setp(ax_dos.get_yticklabels(), visible=False)

plt.suptitle("Quartz Band Structure and DOS")
plt.show()
