import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

plt.rcParams["figure.dpi"] = 150
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["figure.figsize"] = (8, 6)

# #########################
# ---- PHONON DISPERSION ---
# #########################
data = np.loadtxt("si.ph-disp.gp")

k = data[:, 0]
nbands = data.shape[1] - 1

# #########################
# -------- DOS ------------
# #########################
freq, dos, pdos1, pdos2 = np.loadtxt("si.ph-dos.dat", unpack=True)

# #########################
# -------- PLOTTING -------
# #########################
fig = plt.figure()
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05)

ax_band = plt.subplot(gs[0])
ax_dos  = plt.subplot(gs[1], sharey=ax_band)

# ---- Plot phonon branches ----
for band in range(nbands):
    ax_band.plot(k, data[:, band+1], linewidth=1, color='k')

# High symmetry points (ADJUST THESE to your matdyn input)
# Example for L-Γ-X-U-Γ path with 20 points each:
sym_indices = [0, 50, 100, 150, 200,250]

for i in sym_indices:
    ax_band.axvline(k[i], linewidth=0.7, color='k', alpha=0.5)

ax_band.set_xticks([k[i] for i in sym_indices])
ax_band.set_xticklabels(['$\Gamma$','X', 'W','K','$\Gamma$',  'L'])

ax_band.set_ylabel(r"Frequency (cm$^{-1}$)")
ax_band.set_xlim(k.min(), k.max())

# #########################
# ----- plotting DOS ------
# #########################

#ax_dos.plot(dos, freq, color='k', linewidth=1.2)
ax_dos.plot(dos, freq, color='k', lw=1.2, label='Total')
ax_dos.plot(pdos1, freq, '--', lw=1, label='Si (1)')
ax_dos.plot(pdos2, freq, ':', lw=1, label='Si (2)')


# Remove duplicate y tick labels
plt.setp(ax_dos.get_yticklabels(), visible=False)
ax_dos.set_xlabel("DOS")

#plt.suptitle("Silicon Phonon Dispersion and DOS")
#plt.tight_layout()
plt.show()