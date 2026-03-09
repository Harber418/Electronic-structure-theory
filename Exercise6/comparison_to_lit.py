import numpy as np 
import matplotlib.pyplot as plt

c = 2.99792458e10 
plt.figure()
#################### espresso 
data = np.loadtxt("si.ph-disptrue.gp")

nbands = data.shape[1] - 1
for band in range(nbands):
    plt.plot(data[:, 0], data[:, band + 1], linewidth=1, alpha=0.5, color='k')

########################## lit data
for i in range(1,11,1):

    filename = f"si_li{i}.gp"
    lit = np.loadtxt(filename)
    if i <=3:
        x = lit[:, 0] *data[50,0]
        #turn frequency from THz to per cm 
        frequency = lit[:, 1] *1e12/ c
        plt.scatter(x, frequency, marker="o",color="r",s=15)
    elif 8> i >3:
        x = data[100,0] -lit[:, 0]*(data[100,0]-data[50,0])

        #turn frequency from THz to per cm 
        frequency = lit[:, 1] *1e12/ c
        plt.scatter(x, frequency, marker="o",color="r",s=15)
    else:
        x = data[100,0] + 2*lit[:,0]*(data[150,0] - data[100,0])
        #turn frequency from THz to per cm 
        frequency = lit[:, 1] *1e12/ c
        plt.scatter(x, frequency, marker="o",color="r",s=15)
#######################################plotting
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



