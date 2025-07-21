import numpy as np
import matplotlib.pyplot as plt 
import os 
import sys



path    = sys.argv[1]
time,polymeric,monomeric   = np.loadtxt(os.path.join(path,"Log.dat"),usecols=(0,3,4)).T



plt.figure(figsize=(5,5))
plt.plot(time,polymeric)
plt.ylabel("Polymeric Mass",fontsize=15)
plt.xlabel("Time (Sec)",fontsize=15)
plt.tight_layout()
plt.tick_params(labelsize=15)
plt.tight_layout()
plt.show()
