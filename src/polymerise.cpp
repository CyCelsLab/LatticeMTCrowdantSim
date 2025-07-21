#include "class.hpp"        
#include "random.hpp"


void System :: calculateMonoAxis(complexMonomer& i, int newx, int newy, polymer& p){
    int distx,disty,adistx,adisty;        
            
            adisty = abs(i.y - newy);
            adistx = abs(i.x - newx);

            distx  = i.x - newx;
            disty  = i.y - newy;

        if (adistx==0){
            if(adisty>(i.size+1)){
                
                p.h1x = 0;
                p.h2x = 0;
                p.h1y = std :: copysign(1,-disty);
                p.h2y = std :: copysign(1,disty);
            
            }else{

                p.h1x = 0;
                p.h2x = 0;
                p.h1y = std :: copysign(1,disty);
                p.h2y = std :: copysign(1,-disty);
            
            }
        }
        else if(adisty==0){
            if(adistx>(i.size+1)){

                p.h1x = std :: copysign(1,-distx);
                p.h2x = std :: copysign(1,distx);
                p.h1y = 0;
                p.h2y = 0;

            }else{

                p.h1x = std :: copysign(1,distx);
                p.h2x = std :: copysign(1,-distx);
                p.h1y = 0;
                p.h2y = 0;

            }
        }else{

            std :: cout << "Problem hai bhai" << std :: endl;
        }


}


// complex mono-mono fusion 
void System :: complexMonoFusion(complexMonomer& i, complexMonomer& j, int newx, int newy){


    if(i.name==j.name){
        std :: cout << timeInternal << " " << i.name << std :: endl;
        exit(-1);
    }

    // polymer id for new polymer     
    std :: set<int> :: iterator id = IDS.begin(); 

    // new polymer 
    polymer newPolymer;

    // updating monomer i status as it will be incorporated into a polymer 
    i.head = true;
    i.mono = false;
    i.PID  = *id;

    // updating monomer j status as it will be incorporated into a polymer 
    j.head = true;
    j.mono = false;
    j.PID  = *id;


    // updating same information in the lattice 
    fillBlockLattice(j.x,j.y,j.size,0);

    // updating position of monomer j 
    j.x = newx;
    j.y = newy; 

    fillBlockLattice(j.x,j.y,j.size,j.name);
                    
    // filling details in the polymer 
    newPolymer.PID = *id;
    newPolymer.seq = {i.name,j.name};
    newPolymer.length = 2;
    newPolymer.head1 = i.name;
    newPolymer.head2 = j.name;

    // inserting the polymer in the polymer container 
    polymerContainer.insert({*id,newPolymer});
    
    // id is removed from IDs set to prevent duplicate key allocation
    IDS.erase(id);
        
    // update the mass statistics 
    polymerMass+=pow(2*i.size+1,2)*2;
    monomerMass-=pow(2*i.size+1,2)*2;
    //std :: cout << "end" << i.name << " " << j.name << std :: endl;

}


void System :: complexMonoFusion_(complexMonomer& i, complexMonomer& j, int newx, int newy){

    if(i.name==j.name){
        std :: cout << timeInternal << " " << i.name << std :: endl;
        exit(-1);
    }

    // polymer id for new polymer     
    std :: set<int> :: iterator id = IDS.begin(); 

    // new polymer 
    polymer newPolymer;

    // updating monomer i status as it will be incorporated into a polymer 
    i.head = true;
    i.mono = false;
    i.PID  = *id;

    // updating monomer j status as it will be incorporated into a polymer 
    j.head = true;
    j.mono = false;
    j.PID  = *id;



    // updating same information in the lattice 
    fillBlockLattice(j.x,j.y,j.size,0);

    // updating position of monomer j 
    j.x = newx;
    j.y = newy; 

    fillBlockLattice(j.x,j.y,j.size,j.name);

                    
    // filling details in the polymer 
    newPolymer.PID = *id;
    newPolymer.seq = {i.name,j.name};
    newPolymer.length = 2;
    newPolymer.head1 = i.name;
    newPolymer.head2 = j.name;

    // setting polymer-axis
    calculateMonoAxis(i,newx,newy,newPolymer);

    // inserting the polymer in the polymer container 
    polymerContainer.insert({*id,newPolymer});
    
    // id is removed from IDs set to prevent duplicate key allocation
    IDS.erase(id);
        
    // update the mass statistics 
    polymerMass+=pow(2*i.size+1,2)*2;
    monomerMass-=pow(2*i.size+1,2)*2;



}



// function for checking monomers in the vicinity and reacting them. 
void System :: complexMonoFusionCheck(complexMonomer& i){

    // tempvector and tempset will contain the details about neighbours of a monomers
    // which i can possibly bond with it 

    std :: vector<int> tempvector;
    std :: set<int>    tempset; 
    

    int newx,newy;

    // name of the monomers present in the vicinity of monomer i 
    neighbours(i.x,i.y,i.size,tempset);


    for(auto &j : tempset){
        tempvector.push_back(j);
    }

    // shuffle the neighbour to get rid of any bias 
    std :: shuffle(tempvector.begin(),tempvector.end(),engine);

    // loop over all neighbours 
    for(auto &k : tempvector){

       // std :: cout <<11  << std :: endl;

        int x = complexMonomerContainer[k-1].x;
        int y = complexMonomerContainer[k-1].y;

        // get the id of monomer from coordinates 
        int tempMonomerId = k;

        // procced if neighbour is in monomeric form 
        if(complexMonomerContainer[tempMonomerId-1].mono==true){

            // microscopic rate constant 
            if(random<monoBindProb){

                // provide the orientation in which two monomer can bond 
                complexMonoMonoBindingOrientation(i.x, i.y, x, y, i.size, i.size, newx, newy);

                // check if that site is empty 
                if(checkCanMoveBlock(newx,newy,i.size,tempMonomerId)){

                    // polymer formation by monomeric fusion 
                    complexMonoFusion_(i,complexMonomerContainer[tempMonomerId-1], newx,newy);        
                    // loop breaks off 
                    break;
                
                }
            }
        // if monomers present in the visinity is a part of polymer than this block proceeds
        }else if(complexMonomerContainer[tempMonomerId-1].head==true){

            // microscopic rate constant 
            float tempRate  = 0.000000; 
            
            if(polymerContainer[complexMonomerContainer[tempMonomerId-1].PID].length<nucleusSize){
                tempRate = monoBindRate * i.reactionDt;//still not reached threold (nucleating)
            }else{
                tempRate = bindRate * i.reactionDt;//after nucleus formed - elgongating
            }


            if(random<tempRate){
                //std :: cout << "boom" << std :: endl;
                // fusion of monomer and polymer 
                int tempx, tempy, head, PID;
                PID = complexMonomerContainer[tempMonomerId-1].PID;
                complexMonoPolyBindingOrientation(i,complexMonomerContainer[tempMonomerId-1],tempx,tempy, head);

                if(checkCanMoveBlock(tempx,tempy,i.size,i.name)){
                    //std :: cout << "boom2 " << polymerContainer[complexMonomerContainer[tempMonomerId-1].PID].length << std :: endl;
                    complexMonoPolyFusion(i,complexMonomerContainer[tempMonomerId-1],head,PID, tempx,tempy);    
                    break;
                }
                
            }
        }else{

        }
    }
}


void System :: complexPolyFusionCheck(complexMonomer& i){
                

    // tempvector and tempset will contain the details about neighbours of a monomers
    // which i can possibly bond with it 
    std :: vector<int> tempvector;
    std :: set<int>    tempset; 
    
    //int newx,newy;

    // name of the monomers present in the vicinity of monomer i 
    neighbours(i.x,i.y,i.size,tempset);

    for(auto &i : tempset){
        tempvector.push_back(i);
    }


    // shuffle the neighbour to get rid of any bias 
    std :: shuffle(tempvector.begin(),tempvector.end(),engine);

    // loop over all neighbours 
    for(auto &k : tempvector){

        // get the id of monomer from coordinates 
        int tempMonomerId = k;

        // procced if neighbour is in monomeric form 
        if(complexMonomerContainer[tempMonomerId-1].mono==true){

            // microscopic rate constant 
            float tempRate  = 0.000000; 

            if(polymerContainer[i.PID].length<nucleusSize){
                tempRate = monoBindRate * i.reactionDt;
            }else{
                tempRate = bindRate * i.reactionDt;
            }
            if(random<tempRate){
                //std :: cout << "boom" << std :: endl;

                // fusion of monomer and polymer 
                int tempx, tempy, head, PID;
                PID = i.PID;
                complexMonoPolyBindingOrientation_(complexMonomerContainer[tempMonomerId-1],i,tempx,tempy, head);

                if(checkCanMoveBlock(tempx,tempy,i.size,i.name)){
                    //std :: cout << "boom2 " << polymerContainer[i.PID].length << std :: endl;
                    complexMonoPolyFusion(complexMonomerContainer[tempMonomerId-1],i,head,PID, tempx,tempy);    
                    break;
                }
            
            }

        }else if(complexMonomerContainer[tempMonomerId-1].head==true){

                // polymeric polymeric fusion not implemented.

        }else{

        }
    }
}

   
void System ::  complexMonoMonoBindingOrientation(int ix, int iy, int jx, int jy, int  isize, int jsize, int &newx, int &newy){

    /*
    inputs to the function  
    
    ix is the x coordinate of monomer i 
    jx is the x coordinate of monomer j

    iy is the y coordinate of monomer i
    jy is the y coordinate of monomer j
    newx new x coordinate of monomer j 
    newy new y coordinate of monomer j
    */

    float compx,compy;

    vectorComp(ix,iy,jx,jy,compx,compy);
    
    // distance in x and y direction 
    int xdist  = ix - jx;                 
    int ydist  = iy - jy;
    
    // absolute distance in x and y direction 
    int absxdist = abs(xdist);
    int absydist = abs(ydist);


    if(fabs(compx+compy)==1.0){
        newx = ix + (isize + jsize + 1)*compx;
        newy = iy + (isize + jsize + 1)*compy;
        return;
    }

    if(absxdist<absydist){
        if(ydist>=0){
        
            newx = ix; 
            newy = iy - isize - jsize - 1; 
            //std :: cout << 1 << std :: endl;
        
        }else if(ydist<0){
        
            newx = ix;
            newy = iy + isize + jsize + 1;
            //std :: cout << 2 << std :: endl;

        }
    }else if(absxdist>absydist){
        if(xdist>=0){
        
            newx = ix - isize - jsize - 1;
            newy = iy; 
            //std :: cout << 3 << std :: endl;

        
        }else if(xdist<0){
        
            newx = ix + isize + jsize + 1;
            newy = iy; 
            //std :: cout << 4 << std :: endl;

        }
    }else if(absxdist==absydist){

        if(random<.5){
            if(ydist>=0){
            
                newx = ix; 
                newy = iy - isize - jsize - 1; 
                //std :: cout << 5 << std :: endl;

            
            }else if(ydist<0){
            
                newx = ix;
                newy = iy + isize + jsize + 1;
                //std :: cout << 6 << std :: endl;

            }
        }else{
            if(xdist>=0){
            
                newx = ix - isize - jsize - 1;
                newy = iy; 
                //std :: cout << 7 << std :: endl;

            }else if(xdist<0){

                newx = ix + isize + jsize + 1;
                newy = iy; 
                //std :: cout << 8 << std :: endl;

            }
        }
    }
    return;
}

void System ::  complexMonoMonoBindingOrientation_(int ix, int iy, int jx, int jy, int  isize, int jsize, int &newx, int &newy){

    /*
    inputs to the function  
    
    ix is the x coordinate of monomer i 
    jx is the x coordinate of monomer j

    iy is the y coordinate of monomer i
    jy is the y coordinate of monomer j
    newx new x coordinate of monomer j 
    newy new y coordinate of monomer j
    */

    float compx,compy;

    vectorComp(ix,iy,jx,jy,compx,compy);
    
    // distance in x and y direction 
    int xdist  = ix - jx;                 
    int ydist  = iy - jy;
    
    // absolute distance in x and y direction 
    int absxdist = abs(xdist);
    int absydist = abs(ydist);


    if(fabs(compx+compy)==1.0){
        newx = ix + (isize + jsize + 1)*compx;
        newy = iy + (isize + jsize + 1)*compy;
        return;
    }


    if(absxdist<absydist){
        if(ydist>=0){
        
            newx = ix; 
            newy = iy - isize - jsize - 1; 
            //std :: cout << 1 << std :: endl;
        
        }else if(ydist<0){
        
            newx = ix;
            newy = iy + isize + jsize + 1;
            //std :: cout << 2 << std :: endl;

        }
    }else if(absxdist>absydist){
        if(xdist>=0){
        
            newx = ix - isize - jsize - 1;
            newy = iy; 
            //std :: cout << 3 << std :: endl;

        
        }else if(xdist<0){
        
            newx = ix + isize + jsize + 1;
            newy = iy; 
            //std :: cout << 4 << std :: endl;

        }
    }else if(absxdist==absydist){

        if(random<.5){
            if(ydist>=0){
            
                newx = ix; 
                newy = iy - isize - jsize - 1; 
                //std :: cout << 5 << std :: endl;

            
            }else if(ydist<0){
            
                newx = ix;
                newy = iy + isize + jsize + 1;
                //std :: cout << 6 << std :: endl;

            }
        }else{
            if(xdist>=0){
            
                newx = ix - isize - jsize - 1;
                newy = iy; 
                //std :: cout << 7 << std :: endl;

            }else if(xdist<0){

                newx = ix + isize + jsize + 1;
                newy = iy; 
                //std :: cout << 8 << std :: endl;

            }
        }
    }
    return;
}


// only end-on polymerisation
void System :: complexMonoPolyBindingOrientation(complexMonomer& m , complexMonomer& n, int& tempx,int& tempy, int& head){

        // m is the free monomer, while n is a part of a polymer 
        int PID = n.PID;

        // head1 and head2 of the polymer 
        int head1   = polymerContainer[PID].head1;
        int head2   = polymerContainer[PID].head2;
        

        // x and y coordinates of the head1
        int head1x  = complexMonomerContainer[head1-1].x;
        int head1y  = complexMonomerContainer[head1-1].y;
        
        // x and y coordinate of the head2
        int head2x  = complexMonomerContainer[head2-1].x;
        int head2y  = complexMonomerContainer[head2-1].y;

        // varibles for storing components of the unit vector 
        float normx=0,normy=0;

        // either head1 or head1
        whichHead(head1,head2,n.name,head);

        // unit vector of filment 
        vectorComp(head1x,head1y,head2x,head2y,normx,normy);

        if(head==1){
            // one of them will be zero as polymerisation can only happen in basis direction 
            tempx = -(normx*(2*n.size+1))+n.x;
            tempy = -(normy*(2*n.size+1))+n.y;
        }else if(head==2){

            tempx = (normx*(2*n.size+1))+n.x;
            tempy = (normy*(2*n.size+1))+n.y;
        }
        pbcX(tempx);
        pbcY(tempy);
    
}

// only end-on polymerisation
void System :: complexMonoPolyBindingOrientation_(complexMonomer& m , complexMonomer& n, int& tempx,int& tempy, int& head){

        // m is the free monomer, while n is a part of a polymer 
        int PID = n.PID;

        // head1 and head2 of the polymer 
        int head1   = polymerContainer[PID].head1;
        int head2   = polymerContainer[PID].head2;
        
        // either head1 or head1
        whichHead(head1,head2,n.name,head);

        // unit vector of filment 
        //vectorComp(head1x,head1y,head2x,head2y,normx,normy);
        
        if(head==1){
            // one of them will be zero as polymerisation can only happen in basis direction 
            tempx = (polymerContainer[PID].h1x*(2*n.size+1))+n.x;
            tempy = (polymerContainer[PID].h1y*(2*n.size+1))+n.y;
        }else if(head==2){
            tempx = (polymerContainer[PID].h2x*(2*n.size+1))+n.x;
            tempy = (polymerContainer[PID].h2y*(2*n.size+1))+n.y;
        }
     
        pbcX(tempx);
        pbcY(tempy);
    
}
void System ::  complexMonoPolyFusion(complexMonomer& m ,complexMonomer& n,int  head, int PID,int tempx,int tempy){
                

                n.head = false;
                
                m.head = true;
                m.mono = false;
                m.PID  = n.PID;


                fillBlockLattice(m.x,m.y,m.size,0);

                m.x    = tempx;
                m.y    = tempy;

                fillBlockLattice(m.x,m.y,m.size,m.name);


                if(head==1){

                    polymerContainer[PID].seq.push_front(m.name);
                    polymerContainer[PID].head1 = m.name;
                    
                }else if(head==2){

                    polymerContainer[PID].seq.push_back(m.name);
                    polymerContainer[PID].head2 = m.name;
                }
                polymerContainer[PID].length+=1;
        
                polymerMass+= pow(2*m.size+1,2);
                monomerMass-= pow(2*m.size+1,2);
}


void System :: whichHead(int& head1, int& head2, int& name, int& head){

        if(complexMonomerContainer[head1-1].name==name){
        
            head=1;

        }else if(complexMonomerContainer[head2-1].name==name){

            head=2;

        }else{

            std :: cout << complexMonomerContainer[head2-1].name << " "
            << complexMonomerContainer[head1-1].name << " "
            << name << std :: endl;
            std :: cout << "Program failed !!!!" << std :: endl;
            exit(-1);
        }
}

// function to check monomers present in the vicinity of monomers of size, size and present in at site, x,y. 
void System :: neighbours(int x, int y, int size, std :: set<int> &tempset){

    int side = pow(2*size+2,1);
    int tempx = x-size-2;
    int tempy = y+size+1;
    int counter = 0;

    
    int tx = 1;
    int ty = 0;

    int H1,H2,tempH;


    // loop for traversing the perimeter outer to a monnomer 
    for(auto i = 0;i <((side)*4);i++){

        if(i%(side)==0){

            counter+=1;

            switch(counter){

                case 2:
        
                    tempx+=(tx);
                    tempy+=(ty);

                    tx =  0; 
                    ty = -1;
                    break;

                case 3:

                    tempx+=(tx);
                    tempy+=(ty);

                    tx = -1; 
                    ty =  0;
                    break;

                case 4:
                    tempx+=(tx);
                    tempy+=(ty);

                    tx =  0; 
                    ty =  1;
                    break;

                default:

                    tempx+=(tx);
                    tempy+=(ty);
                break;

            }
        }else{   

            tempx+=(tx);
            tempy+=(ty);
        }
        
        // PBC conditions 
        pbcX(tempx);
        pbcY(tempy);
        
        // Only valid monomers will be put in tempset. As number of monomers is from 1 to 
        // total monomers. 


        if((tempx>=0)&&(tempy>=0)&&(tempx<NX)&&(tempy<NY)){
            if((lattice[tempx][tempy]!=0)&&(lattice[tempx][tempy]<=nComplexMonomers)){

                
                int tempid   = lattice[tempx][tempy];

                // either neighbour should be in monomeric state 
                if(complexMonomerContainer[tempid-1].mono){                
                    tempset.insert(tempid);
                
                // or it a head of a monomer
                }else if(complexMonomerContainer[tempid-1].head){

                    if(getDpolType()==1){

                        H1 = polymerContainer[complexMonomerContainer[tempid-1].PID].head1;
                        H2 = polymerContainer[complexMonomerContainer[tempid-1].PID].head2;
                        
                        whichHead(H1,H2,complexMonomerContainer[tempid-1].name,tempH);
                        
                        if(tempH==2){
                            tempset.insert(tempid);
                        }

                    }else{

                        tempset.insert(tempid);

                    }

                }

            }
        }
    }
}
