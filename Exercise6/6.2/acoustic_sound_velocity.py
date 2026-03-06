import numpy as np
import matplotlib.pyplot as plt

#####################################################
#------looking for the sound velocity of silicon 
#################################################

#we need to find the slope of the acoustic branches near gamma 

#read in the data
data = np.loadtxt("si.ph-disp7.gp")
#the first 3 bands are all we are interested in 
nbands = data.shape[1] - 1
names = ["transverse","transverse 2","longitudinal"]
colors = ['r', 'g', 'b']
for band in range(1,4,1):
    plt.plot(data[:, 0], data[:, band ], linewidth=1, alpha=0.5, color='k')
    #we want to extract the first 10 points of the acoustic branches to find the gradiant
    x = data[:2, 0]
    y = data[:2, band ]
    coeffs = np.polyfit(x, y, 1)
    slope = coeffs[0]
    if band == 2:
        pass
    else:
        plt.plot(data[:20, 0], np.polyval(coeffs, data[:20, 0]), linestyle='--', color=colors[band-1],label=f"{names[band-1]} slope = {slope:.2f}",alpha = 0.5)
    #find the sound velocity 
    #y data is in cm^-1 this is Kayser
    #x data is in reduced coordinates, and needs to be converted to meters 
    #
    v = slope * 1e4 * 1e-10 / (1.0545718e-34 / 1.66053906660e-27) # convert from cm/s to m/s and from eV to J
    #v = slope * 2* np.pi * 2.99*10**10
    print(f"{names[band-1]} sound velocity: {v:.2f} m/s")
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
plt.legend()
plt.show()
#  0.000  0.000  0.000  50
 # 0.500  0.000  0.500  50
 # 0.500  0.250  0.750  50
 # 0.375  0.375  0.750  50
  #0.000  0.000  0.000  50
  #0.500  0.500  0.500  50