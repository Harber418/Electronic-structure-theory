import numpy as np
import matplotlib.pyplot as plt

# Load phonon DOS

freq, dos, pdos1, pdos2 = np.loadtxt("si.ph-dos.dat", unpack=True)

plt.figure(figsize=(6,6))

plt.plot(freq, dos, color='k', lw=1.2, label='Total')
plt.plot(freq, pdos1, '--', lw=1, label='Si (1)')
plt.plot(freq, pdos2, ':', lw=1, label='Si (2)')

plt.xlabel(r'Frequency (cm$^{-1}$)')
plt.ylabel(r'Phonon DOS (states/cm$^{-1}$/u.c.)')

plt.xlim(0, freq.max())
plt.ylim(bottom=0)

plt.legend(frameon=False)
plt.tight_layout()
plt.show()