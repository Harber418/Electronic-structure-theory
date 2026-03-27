import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt

def read_data():
    #o--o   anlge o--h volume
    ice1= np.array([2.73726,107.5909,1.05445,1736.0838])
    ice1_down = [2.54018,107.5909,0.97853,1387.4413]
    ice1_up = [2.92613,107.5909,1.12721,2120.8200]
    ice2 = [2.76551,99.8624,1.01312,2053.1832]
    ice2_down =[2.56640,99.8624,0.94017,1640.8604]
    ice2_up =[ 2.95633,99.8624,1.08302,2508.1922]
    ice8 =[2.89541,105.5976,0.96857,1008.5161]
    ice8_down =[2.68694,105.5976,0.89883,805.9846]
    ice8_up =[3.09519,105.5976,1.03540,1232.0148]
    #van 
    ice1_van =[2.75760,106.6896,0.99815]
    ice1_van_up =[2.94787,06.6896,1.06702]
    ice1_van_down =[2.55905,106.6896,0.92628]
    ice2_van =[2.77777,106.8328,0.99188]
    ice2_van_up =[2.96944,106.8328,1.06032]
    ice2_van_down =[2.57777,106.8328,0.92047]
    ice8_van =[2.89540,105.9522,0.99017]
    ice8_van_down =[2.68693,105.9522,0.91888]
    ice8_van_up =[3.09518,105.9522,1.05849]

    #plot the 0--0 bond length 
    plt.figure()
    volume = [ice1_down[3],ice1[3],ice1_up[3]]*0.529177 ** 3
    length = [ice1_down[0],ice1[0],ice1_up[0]]
    volume2 = [ice2_down[3],ice2[3],ice2_up[3]]*0.529177 ** 3
    length2 = [ice2_down[0],ice2[0],ice2_up[0]]
    volume8 = [ice8_down[3],ice8[3],ice8_up[3]]*0.529177 ** 3
    length8 = [ice8_down[0],ice8[0],ice8_up[0]]

    plt.plot(volume,length,color = "red",label="ice1")
    plt.plot(volume2,length2,color = "orange",label="ice2")
    plt.plot(volume8,length8,color = "gold",label="ice8")
    plt.legend()
    plt.show()

    #plot o--h
    plt.figure()
    volume = [ice1_down[3],ice1[3],ice1_up[3]]*0.529177 ** 3
    length = [ice1_down[2],ice1[2],ice1_up[2]]
    volume2 = [ice2_down[3],ice2[3],ice2_up[3]]*0.529177 ** 3
    length2 = [ice2_down[2],ice2[2],ice2_up[2]]
    volume8 = [ice8_down[3],ice8[3],ice8_up[3]]*0.529177 ** 3
    length8 = [ice8_down[2],ice8[2],ice8_up[2]]

    plt.plot(volume,length,color = "red",label="ice1")
    plt.plot(volume2,length2,color = "orange",label="ice2")
    plt.plot(volume8,length8,color = "gold",label="ice8")
    plt.legend()
    plt.show()


    
if __name__ == "__main__":
    read_data()