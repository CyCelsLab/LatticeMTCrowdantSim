#include "class.hpp"
#include "random.hpp"


void System :: preseeds(int iTotalPreSeeds, int length, int size){

    System :: preseededPolymerSize = length;
    
    int counter = 0;
    while(counter<iTotalPreSeeds){

        int count = 0, tempx = random*NX,tempy = random*NY,normx=0,normy=0, chance = 4*random;


         switch(chance){
            case 0:
                normx=1*(2*size+1);
                normy=0;
                break;
            case 1:
                normx=-1*(2*size+1);
                normy=0;
                break;
            case 2:
                normy=1*(2*size+1);
                normx=0;
                break;
            case 3: 
                normy=-1*(2*size+1);
                normx=0;
                break;
        }

        for(auto i = 0;i<length;i++){

            int tempxx=tempx+normx*i;
            pbcX(tempxx);           
            int tempyy=tempy+normy*i;
            pbcY(tempyy);  

            //std :: cout << size << std :: endl;
            if(checkEmptyBlock(tempxx,tempyy,size)){
                //if(checkBlockInLattice(tempx,tempy,size)){                
                    count+=1;
                //}
            }
        }


        


        if(count==length){

            int tempcount=0;
            for(int i =0;i<NX;i++){
                for(int j=0;j<NY;j++){
                    if(lattice[i][j]){
                    tempcount+=1;
                    }
                }
            }


           // std :: cout << counter << " " << count << " " << tempcount<< std :: endl;

            polymer tempPolymer;
            
            tempPolymer.head1 = nComplexMonomers + 1;
            tempPolymer.length = length;
            tempPolymer.PID    = counter+1;

            for(auto i=0;i<length;i++){
            
                complexMonomer tempMonomer;
                tempMonomer.name = nComplexMonomers+1;
                nComplexMonomers++;

                int tempxx=tempx+normx*i;
                pbcX(tempxx);           
                int tempyy=tempy+normy*i;
                pbcY(tempyy);  


                tempMonomer.x = tempxx;
                tempMonomer.y = tempyy;
                tempMonomer.mono=false;

                tempMonomer.PID=(counter+1);
                tempMonomer.head=false;
                tempMonomer.size= size;  
                
                fillBlockLattice(tempMonomer.x,tempMonomer.y,tempMonomer.size,tempMonomer.name);                
                
                if(i==0){tempMonomer.head=true;}
                if(i==(length-1)){tempMonomer.head=true;}
                
                complexMonomerContainer.push_back(tempMonomer);
                tempPolymer.seq.push_back(tempMonomer.name);
                totalMass+=pow(2*size+1,2);
                //monomerMass+=pow(2*size+1,2);
            
                //std :: cout << 3 << " " << i  << std :: endl;

                randomEntity.push_back(tempMonomer.name);

            
            }
            
            tempPolymer.head2 = nComplexMonomers;
            
            if(normx==0){
                tempPolymer.h1x   = 0;
                tempPolymer.h2x   = 0;
            }else{
                tempPolymer.h1x    = std :: copysign(1,-normx);
                tempPolymer.h2x    = std :: copysign(1, normx);
            }

            if(normy==0){
                tempPolymer.h1y   = 0;
                tempPolymer.h2y   = 0;
            }else{
                tempPolymer.h1y    = std :: copysign(1,-normy);
                tempPolymer.h2y    = std :: copysign(1, normy);
            }

            polymerContainer.insert({tempPolymer.PID,tempPolymer});
            
            polymerMass+=pow(2*size+1,2)*length;
            
            counter++;

            preseededPolymers++;


            tempcount=0;
            for(int i =0;i<NX;i++){
                for(int j=0;j<NY;j++){
                    if(lattice[i][j]){
                    tempcount+=1;
                    }
                }
            }
            //std :: cout << counter << " " << count << " " << tempcount<<  " " << counter*225 << std :: endl;


        }
        
    }
}
