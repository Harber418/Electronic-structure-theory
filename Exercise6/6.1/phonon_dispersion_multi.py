import numpy as np
import matplotlib.pyplot as plt
for i in range(1,11,1):
    data = np.loadtxt(f"si.ph-disp{i}.gp")

    nbands = data.shape[1] - 1
    colours = ["red", "blue", "green", "orange", "purple", "cyan", "magenta", "brown", "pink", "gray"]
    #plt.plot(data[:, 0], data[:, 1], linewidth=1, alpha=0.5, color=colours[band])
    for band in range(nbands):
        plt.plot(data[:, 0], data[:, band + 1], linewidth=1, alpha=0.5, color=colours[band])
        plt.axvline(x=data[0, 0], linewidth=0.5, color=colours[band], alpha=0.01,label=f"q = {i}")
        
    # High symmetry k-points (check matdyn.GaAs.in)
plt.axvline(x=data[0, 0], linewidth=0.5, color='k', alpha=0.5)
plt.axvline(x=data[50, 0], linewidth=0.5, color='k', alpha=0.5)
plt.axvline(x=data[100, 0], linewidth=0.5, color='k', alpha=0.5)
plt.axvline(x=data[150, 0], linewidth=0.5, color='k', alpha=0.5)
plt.axvline(x=data[200, 0], linewidth=0.5, color='k', alpha=0.5)
plt.axvline(x=data[250, 0], linewidth=0.5, color='k', alpha=0.5)
plt.xticks(ticks= [0, data[50, 0], data[100, 0], data[150, 0],data[200,0], data[-1, 0]], \
           labels=['$\Gamma$','X', 'W','K','$\Gamma$',  'L'])
plt.ylabel("Frequency (cm$^{-1}$)")
plt.xlim(data[0, 0], data[-1, 0])
plt.ylim(0, )
plt.show()
