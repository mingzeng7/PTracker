import numpy as np
import matplotlib.pyplot as plt
import h5py

with h5py.File('output.h5','r') as h_file:
    data = h_file["0"]
    #plt.plot(data[:,0], data[:,1], label='1')
    plt.plot(data[:,0], data[:,2], label='2')
    #plt.plot(data[:,0], data[:,3], label='3')
    #plt.plot(data[:,0], data[:,4], label='4')
    plt.legend()
    plt.show()
