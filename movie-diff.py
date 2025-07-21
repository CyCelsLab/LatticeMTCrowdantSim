import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import matplotlib


# input file directory should be given via command line argument 
src = sys.argv[1]


cmap = matplotlib.colors.ListedColormap(["white","green","red","blue"], name='from_list', N=None)

#printing iteration 
pritingfile = 500
#ermination criterion 
endpoint    = 2001

with open(os.path.join(src,"Log.dat"),"r") as fileloop:
    k = 0 
    for line in fileloop:
        k+=1

        if k==1:
            sizem = int(line.split()[3])
                    
        if k==2:
            sizeX = int(line.split()[-1])

        if k==3:
            sizeY = int(line.split()[-1])
         
        if k==16:
            sizec = int(line.split()[4])
         
         

        if(line.startswith("#")==0): break
        
print(sizeX,sizeY)



with open(os.path.join(src,"Monomer.dat"),"r") as fileloop:
    k= -1
    mdata = []
    for line in fileloop:
        if line[:6] == "#start":

            line = next(fileloop) 
            k+=1
            if k%pritingfile==0:
                read = 1
                mdata.append([])
            if k == endpoint : break    
           

        if line[0:4] == "#end":
            read = 0

        if read:
            temp = [i for i in map(int,line.split())]
            mdata[-1].append([temp[0],temp[1],temp[2],temp[-1]])

mdata = np.array(mdata[:-1])

with open(os.path.join(src,"Crowdant.dat"),"r") as fileloop:
    k= -1
    cdata = []
    for line in fileloop:
        if line[:6] == "#start":
            k+=1
            line = next(fileloop) 

            if k%pritingfile==0:
                read = 1
                cdata.append([])
            if k == endpoint : break    


        if line[0:4] == "#end":
            read = 0

        if read:
            temp = [i for i in map(int,line.split())]
            cdata[-1].append([temp[0],temp[1],temp[2]])
        
cdata = np.array(cdata[:-1])

for kk in range(mdata.shape[0]):
    tempM = np.zeros((sizeX,sizeY))
    
    size = sizem
    for k in mdata[kk]:

        index,x,y,p = k

        if p > 0 : 
            for i in np.arange(-1,2):
                for j in np.arange(-1,2):
                    tempx = (x+i)%sizeX
                    tempy = (j+y)%sizeX

                    if tempM[tempx,tempy]==0:
                        tempM[tempx,tempy] =  3
                    else:
                        print(k,kk,index,tempM[tempx,tempy],p)
                        break
        else:
            for i in np.arange(-size,size+1):
                for j in np.arange(-size,size+1):
                    tempx = (x+i)%sizeX
                    tempy = (j+y)%sizeY

                    if tempM[tempx,tempy]==0:
                        tempM[tempx,tempy] = 1
                    else:
                        break
    size = sizec
    for k in cdata[kk]:

        index,x,y = k
        for i in np.arange(-size,size+1):
            for j in np.arange(-size,size+1):
                tempx = (x+i)%sizeX
                tempy = (j+y)%sizeY

                if tempM[tempx,tempy]==0:
                    tempM[tempx,tempy] = 2
                    pass
                else:
                    break

    plt.figure(figsize=(10,10))
    plt.pcolormesh(tempM,cmap=cmap)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(src,f"frame_{kk}.png"),format="png")
    plt.close()
