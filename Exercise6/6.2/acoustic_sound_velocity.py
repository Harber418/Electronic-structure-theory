import numpy as np
import matplotlib.pyplot as plt

#####################################################
#------looking for the sound velocity of silicon 
#################################################

#we need to find the slope of the acoustic branches near gamma 

#read in the data
data = np.loadtxt("si.ph-disp.gp")
#the first 3 bands are all we are interested in 
nbands = data.shape[1] - 1
for band in range(1,4,1):
    plt.plot(data[:, 0], data[:, band + 1], linewidth=1, alpha=0.5, color='k')
# High symmetry k-points (check matdyn.gaas.in)

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
#  0.000  0.000  0.000  50
 # 0.500  0.000  0.500  50
 # 0.500  0.250  0.750  50
 # 0.375  0.375  0.750  50
  #0.000  0.000  0.000  50
  #0.500  0.500  0.500  50