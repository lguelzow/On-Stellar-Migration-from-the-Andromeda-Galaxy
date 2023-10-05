import numpy as np
from latqcdtools.base.cleanData import *

# read in goes as speed, lower error, upper error

print("\nLoading data...\n")
data = np.loadtxt("Marchetti_v_data.txt", unpack=True)

vcut = 637.43

print("Clipping data...\n")
data = clipRange(data, col=0, minVal=vcut)

outFile = open('dataMarchetti_clipped_637.d', 'w')

for i in range(len(data[0])):
    outFile.write(str(data[0, i])+"  "+str(data[1, i]) +
                  "  "+str(data[2, i])+"\n")

outFile.close()
