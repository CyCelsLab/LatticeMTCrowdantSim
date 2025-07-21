import numpy as np
import matplotlib.pyplot as plt 
import os 
import sys



allFrame = []
time = []    

with open(os.path.join(sys.argv[1],"Polymer.dat")) as fileloop:    
    for line in fileloop:
        
        if line[:6] == "#start":
            #print(line)
            time.append(float(line.split()[-1]))
            read  = 1
            temp = np.zeros(1000)
            line = next(fileloop)
        if line[:4] == "#end":
            read = 0
            allFrame.append(temp)
        if read:
            tempi = int(line.split()[2])
            temp[tempi-1]+=1
allFrame = np.array(allFrame)


for i in range(1000):
	if(allFrame[:,i].sum()): plt.plot(time,allFrame[:,i],label=f"{i}")
        
plt.xlabel("Time (sec)",fontsize=15)
plt.ylabel("Polymer Number",fontsize=15)
plt.tick_params(labelsize=15)
plt.legend(title="Polymer Length")
plt.tight_layout()
plt.show()
