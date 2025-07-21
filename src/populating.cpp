//populating.cpp
#include"class.hpp"
#include "random.hpp"

//
bool System :: checkBlockInLattice(int& x, int& y, int& size){

        // assuming size of the block to be square
    
        // assuming that the block is in the lattice
        bool flag=1;
        
        // check for x axis lower bound
        if((x-size)<0){
            flag=0;
            return(flag);
        }

        // check for x axis upper bound
        if((x+size)>=NX){
            flag=0;
            return(flag);
        }

        // check for x axis upper bound
        if((y-size)<0){
            flag=0;
            return(flag);
        }

        // check for x axis upper bound
        if((y+size)>=NY){
            flag=0;
            return(flag);
        }

        return(flag);
}

// function for filling the lattice with a specific value
void System :: fillBlockLattice(int& x, int& y, int& size, int value){

        int i=-1,j=-1;

        for(int m=-size;m<=size;++m){
            for(int n=-size;n<=size;++n){

                i = (NX+(x+(m))%NX)%NX;
                j = (NY+(y+(n))%NY)%NY;

                //std :: cout << x << " " << i << " " << y << " " << j << std :: endl;

                lattice[i][j] = value;
            
            }
        }
}

// function to check if block of size x size is empty
bool System :: checkCanMoveBlock(int& x, int& y, int& size, int name){

        // assuming lattice sites are vacant
        bool flag = 1;
        
        int i=0,j=0;

        //std :: cout << x << " " << y << " " << name << std :: endl; 
        for(int m=-size;m<=size;++m){
            for(int n=-size;n<=size;++n){


                i = (NX+(x+(m))%NX)%NX;
                j = (NY+(y+(n))%NY)%NY;

           // std :: cout << i << " " << j << " " << lattice[i][j]<< std :: endl; 

                if(lattice[i][j]>0){

                    if(lattice[i][j]!=name){
                        flag=0; 
                        return(flag);
                    }

                }
            }
        }

       // std :: cout << flag << std :: endl;
        // if all the sites are unoccupied reture true.
        if(flag==1){
            return(flag);
        }else{

            std :: cout << "error is fired" << std :: endl;
            return(0);
        }
}

                               
void System :: populateLattice(int iTotalComplexCrowdants, int csize, int cdiffStepLength, int iTotalComplexMonomers, int msize, int mdiffStepLength, int ipreseeds, int ipreseedSize, int iDistributionPattern){


    int total     = iTotalComplexMonomers + iTotalComplexCrowdants;
    int temp1     = iTotalComplexMonomers, temp2 = iTotalComplexCrowdants+iTotalComplexMonomers, counter = 0; 
    
    int chance    = 0;
    int tempmono  = 0;
    int tempcrow  = 0;
    int totalMonomers = iTotalComplexMonomers + ipreseeds*ipreseedSize;

    monomerStepLength   = mdiffStepLength;
    crowdantStepLength  = cdiffStepLength;



    distributionPattern = iDistributionPattern; // 0 == random // 1 == monomer First // 2 == crowdant First // 

    //System :: preseeds(ipreseeds,ipreSeedSize,mdiff,msize);


    for(auto i=0; i<total;i++){
        
        if(distributionPattern==0){
            
            chance = random*(total-counter);

            if(chance<temp1){
                
                tempmono+=1;
                
                ++nComplexMonomers;
                ++totalEntity;

                putComplexMonomers(nComplexMonomers,msize,mdiffStepLength);

                temp1-=1;
                temp2-=1;
            
            }else if(chance>=temp1&&chance<temp2){
                
                tempcrow+=1;
                ++nComplexCrowdants;
                ++totalEntity;

                putComplexCrowdants(nComplexCrowdants+totalMonomers,csize,cdiffStepLength);    
                temp2-=1;

            }else{
            
                std :: cout << "Allocation Failed" << std :: endl;
                exit(-1);
            
            }       
            counter+=1;
        }else if(distributionPattern==1){


            if(i < temp1){
                
                ++nComplexMonomers;
                ++totalEntity;
                putComplexMonomers(nComplexMonomers,msize,mdiffStepLength);
            
            }else{
                
                ++nComplexCrowdants;
                ++totalEntity;
                putComplexCrowdants(nComplexCrowdants+totalMonomers,csize,cdiffStepLength);
            
            }
 
        }else if(distributionPattern==2){
            

            temp2 = iTotalComplexCrowdants;
            temp1 = iTotalComplexCrowdants + iTotalComplexMonomers;

            if(i < temp2){
            
                ++nComplexCrowdants;
                ++totalEntity;
                putComplexCrowdants(nComplexCrowdants+totalMonomers,csize,cdiffStepLength);
            
            }else{
            
                ++nComplexMonomers;
                ++totalEntity;
                putComplexMonomers(nComplexMonomers,msize,mdiffStepLength);
            
            }
        }
    }     
   // std :: cout << nComplexCrowdants << " " << nComplexMonomers << std :: endl;
}

// function to check if block of size x size is empty
bool System :: checkEmptyBlock(int& x, int& y, int& size){

        // assuming lattice sites are vacant
        bool flag = 1;
        
        int i=0,j=0;

        for(int m=-size;m<=size;++m){
            for(int n=-size;n<=size;++n){


                i = (NX+(x+(m))%NX)%NX;
                j = (NY+(y+(n))%NY)%NY;

                if(lattice[i][j]>0){

                        flag=0; 
                        return(flag);
                }
            }
        }
        return(flag);
}


// taking complex crowdant in the system 

void System :: putComplexCrowdants(int icomplexCrowdantName, int size, int stepLength){

        // size of the complex crowdant will be in terms of moore neighbo
        int tempi,tempj, counter=0;

       // in SI units 
        double radius       = (2*size+1)*gridSpacing*1e-7/2; 
        // in SI units 
        crowdiff    = (4.1 * 1e-21)/(6*M_PI*1e-3*radius); 

        // converting back to length units of 1e-7;
        crowdiff =  crowdiff*1e14;

        // diffusion coefficient from diffusion rate 
        crowdantDiffusionRate = (4*crowdiff)/(stepLength*stepLength * float(gridSpacing* gridSpacing));
        //crowdantDiffusionRate = (4*diff)/(stepLength*stepLength);


        while(counter<1){

            tempi = size + (NX-(2*size))*random;
            tempj = size + (NY-(2*size))*random;


           // tempi = NX/4 + size + (NX/2-(2*size))*random;
           // tempj = NX/4 + size + (NY/2-(2*size))*random;


            if(checkEmptyBlock(tempi,tempj,size)){

                
                complexCrowdant temp;
                
                temp.x = tempi;
                temp.y = tempj;
                temp.diff = crowdiff;
                temp.size = size;
                temp.mass = pow(2*size+1,2);
                temp.waitingTime = -log(random)/crowdantDiffusionRate;

                crowdantMass+=pow(2*size+1,2);
                totalMass+= pow(2*size+1,2);  

                temp.name   = icomplexCrowdantName;

                fillBlockLattice(tempi,tempj,size,temp.name);

                complexCrowdantContainer.push_back(temp);
                counter+=1;
                randomEntity.push_back(temp.name);
               // std :: cout << counter << std :: endl;
            }

        
        }
        
}


void System :: putComplexMonomers(int iComplexMonomerName, int size, int stepLength){

        // std :: cout << iComplexMonomerName << std :: endl;
        // size of the complex crowdant will be in terms of moore neighbour

        int tempi,tempj, counter=0;



        // in SI units 
        double radius       = (2*size+1)*gridSpacing*1e-7/2; 
        // in SI units 
        monodiff    = (4.1 * 1e-21)/(6*M_PI*1e-3*radius); 


        // converting back to length units of 1e-7;
        monodiff =  monodiff*1e14;


        //monomerDiffusionRate = (4*diff)/(stepLength*stepLength);
        //monomerDiffusionRate = (4*diff)/(stepLength*stepLength * float(gridSpacing* gridSpacing));

        monomerDiffusionRate = (4*monodiff)/(stepLength*stepLength * float(gridSpacing* gridSpacing));


        while(counter<1){

            tempi = size + (NX-(2*size))*random;
            tempj = size + (NY-(2*size))*random;

            // tempi = NX/4 + size + (NX/2-(2*size))*random;
            // tempj = NX/4 + size + (NY/2-(2*size))*random;


            if(checkEmptyBlock(tempi,tempj,size)){

                complexMonomer temp;                 

                temp.x = tempi;
                temp.y = tempj;

                temp.size = size;

                temp.mono = true;
                temp.head = false;
                temp.PID  = -1;
                temp.waitingTime = -log(random)/monomerDiffusionRate;
                temp.reactionDt  = temp.waitingTime;

                temp.name   = iComplexMonomerName;
                ++counter;
                fillBlockLattice(tempi,tempj,size,temp.name);
                
                complexMonomerContainer.push_back(temp);
                
                totalMass+= pow(2*size+1,2);  
                monomerMass+=pow(2*size+1,2);
                                
                randomEntity.push_back(temp.name);
               
            }
        }
}

// setting the possible polymer ids 
void System :: setIDS(){
    for(int i=0;i<((nComplexMonomers+1)/2);i++){ 
        IDS.insert(i+1);
    }
}