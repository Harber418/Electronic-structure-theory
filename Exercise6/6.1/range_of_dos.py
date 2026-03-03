import numpy as np
import matplotlib.pyplot as plt
# Load phonon DOS
plt.figure(figsize=(6,6))
for i in range(1, 11, 1 ):
    
    freq, dos, pdos1, pdos2 = np.loadtxt(f"si_dos_{i}.dat", unpack=True)
    dos = dos + i*0.1
    plt.plot(freq, dos, lw=1.2, label=f"q = {i}")
    #plt.plot(freq, pdos1, '--', lw=1, alpha= 0.3)
    #plt.plot(freq, pdos2, ':', lw=1,alpha =0.3)

plt.xlabel(r'Frequency (cm$^{-1}$)')
plt.ylabel(r'Phonon DOS (states/cm$^{-1}$/u.c.)')

plt.xlim(0, freq.max())
plt.ylim(bottom=0)

plt.legend(loc="upper left")
#plt.legend(frameon=False)
plt.tight_layout()
plt.show()