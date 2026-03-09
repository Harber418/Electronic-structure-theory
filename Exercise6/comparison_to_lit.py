import numpy as np 
import matplotlib.pyplot as plt

c = 2.99792458e10 
plt.figure()
for i in range(1,12,1):

    filename = f"si.li{i}.gp"
    data = np.loadtxt(filename)
    x = data[:, 0]
    #turn frequency from THz to per cm 
    frequency = data[:, 1] *1e-12/ c
    plt.plot(x, frequency)



data = np.loadtxt("si.ph-disp.gp")

nbands = data.shape[1] - 1
for band in range(nbands):
    plt.plot(data[:, 0], data[:, band + 1], linewidth=1, alpha=0.5, color='k')
# High symmetry k-points (check matdyn.GaAs.in)
plt.axvline(x=data[0, 0], linewidth=0.5, color='k', alpha=0.5)
plt.axvline(x=data[50, 0], linewidth=0.5, color='k', alpha=0.5)
plt.axvline(x=data[100, 0], linewidth=0.5, color='k', alpha=0.5)
plt.axvline(x=data[150, 0], linewidth=0.5, color='k', alpha=0.5)

plt.xticks(ticks= [0, data[50, 0], data[100, 0], data[-1, 0]], \
           labels=['$\Gamma$','X','$\Gamma$',  'L'])
plt.ylabel("Frequency (cm$^{-1}$)")
plt.xlim(data[0, 0], data[-1, 0])
plt.ylim(0, )
plt.show()
#path for si in fcc 
#  0.000  0.000  0.000  50 gamma
 # 0.500  0.000  0.500  50 x
 # 0.500  0.250  0.750  50 w
 # 0.375  0.375  0.750  50 k
  #0.000  0.000  0.000  50 gamma 
  #0.500  0.500  0.500  50 L



