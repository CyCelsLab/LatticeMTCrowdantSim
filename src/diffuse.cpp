// diffuse.cpp
#include <chrono>
#include "class.hpp"
#include "random.hpp"

// function of diffusing complex crowdants
void System ::complexCrowdantDiffuse(complexCrowdant &i){
    int newX, newY,move;

    // choosing a move by random 
    nextMove(i.x, i.y, crowdantStepLength, newX, newY, move);
    crowTotalCounter[move]+=1;

    // check if move is valid or not
    // if valid then crowdant is translated to new location 
    if (checkCanMoveBlock(newX, newY, i.size, i.name)){
        crowHopCounter[move]+=1;
        int temp = i.name;

        fillBlockLattice(i.x, i.y, i.size, 0);
        fillBlockLattice(newX, newY, i.size, temp);
        if(newX<0||newY<0){
            std :: cout << "code blast crowdant" << std :: endl;
            std :: cout << i.name << " " <<  i.x << " " << i.y << std :: endl;
            std :: cout << i.name << " " <<  newX << " " << newY << std :: endl;
        }
        i.x = newX;
        i.y = newY;
    }

}

void System ::complexMonomerDiffuse(complexMonomer &i){

    int newX, newY, move;
    
    // choosing a move by random 
    nextMove(i.x, i.y, monomerStepLength, newX, newY,move);
    monoTotalCounter[move]+=1;

    // check if move is valid or not
    // if valid then monomer is translated to new location
    if (checkCanMoveBlock(newX, newY, i.size, i.name)){
        monoHopCounter[move]+=1;

        int temp = i.name;

        fillBlockLattice(i.x, i.y, i.size, 0);
        fillBlockLattice(newX, newY, i.size, temp);

        if(newX<0||newY<0){
            std :: cout << "code blast monomer" << std :: endl;
            std :: cout << i.name << " " <<  i.x << " " << i.y << std :: endl;
            std :: cout << i.name << " " <<  newX << " " << newY << std :: endl;
            exit(-1);

        }
        i.x = newX;
        i.y = newY;
    }
}

void System :: equillibrate(int iter){

    equllibriationTime = iter;

    for(auto j = 0 ; j < iter ; j++){

        // randomisation of monomers before iteration 
        std :: shuffle(randomEntity.begin(),randomEntity.end(),engine);

        for (auto &i : randomEntity){


            if(i<=nComplexMonomers){
                // monomer will diffuse iff it is in free state 
                if(complexMonomerContainer[i-1].mono){
                    
                    while(complexMonomerContainer[i-1].waitingTime<dt){

                        complexMonomerDiffuse(complexMonomerContainer[i-1]);
                        complexMonomerContainer[i-1].waitingTime+=-log(random)/monomerDiffusionRate;
                    }

                     complexMonomerContainer[i-1].waitingTime-=dt;




                }
            }else{

                    while(complexCrowdantContainer[i-nComplexMonomers-1].waitingTime<dt){

                        complexCrowdantDiffuse(complexCrowdantContainer[i-nComplexMonomers-1]);
                        complexCrowdantContainer[i-nComplexMonomers-1].waitingTime+=-log(random)/crowdantDiffusionRate;
                    }
                     complexCrowdantContainer[i-nComplexMonomers-1].waitingTime-=dt;


            }
        }
    }
}


// function for evoking diffusion in system
void System ::diffuseandReact(){

    int H1;
    int H2;
    int tempH;



    // randomisation of monomers before iteration 
    std :: shuffle(randomEntity.begin(),randomEntity.end(),engine);

    // looping over all radomised monomers list 
    for (auto &i : randomEntity){


        if(i<=nComplexMonomers){

            // monomer will diffuse iff it is in free state 
            if(complexMonomerContainer[i-1].mono){
                

                    while(complexMonomerContainer[i-1].waitingTime<dt){

                        if (complexMonomerContainer[i-1].mono){

                            complexMonomerDiffuse(complexMonomerContainer[i-1]);

                            complexMonoFusionCheck(complexMonomerContainer[i-1]);

                            reactionDt = -log(random)/monomerDiffusionRate;                            
                            
                            complexMonomerContainer[i-1].reactionDt = reactionDt;
                            complexMonomerContainer[i-1].waitingTime+=reactionDt;                    


                        }else{
                            
                            complexMonomerContainer[i-1].reactionDt = fabs(fmod(complexMonomerContainer[i-1].reactionDt,dt));
                            

                            if(getDpolType()==1){
                        
                                H1 = polymerContainer[complexMonomerContainer[i-1].PID].head1;
                                H2 = polymerContainer[complexMonomerContainer[i-1].PID].head2;
                        
                                whichHead(H1,H2,complexMonomerContainer[i-1].name,tempH);
                                    
                                if(tempH==2)  complexPolyFusionCheck(complexMonomerContainer[i-1]);
                            
                            }else{
                            
                                complexPolyFusionCheck(complexMonomerContainer[i-1]);
                            
                            }

                            complexMonomerContainer[i-1].waitingTime = 0;
                            //complexMonomerContainer[i-1].waitingTime+=-log(random)/monomerDiffusionRate;                    
                            break;
                        
                        }

                    }

                     complexMonomerContainer[i-1].waitingTime-=dt;
    
            }else if(complexMonomerContainer[i-1].head){
            
                    complexMonomerContainer[i-1].reactionDt = dt;
                    
                    if(getDpolType()==1){
                    
                        H1 = polymerContainer[complexMonomerContainer[i-1].PID].head1;
                        H2 = polymerContainer[complexMonomerContainer[i-1].PID].head2;
                    
                        whichHead(H1,H2,complexMonomerContainer[i-1].name,tempH);
                        
                        if(tempH==2) complexPolyFusionCheck(complexMonomerContainer[i-1]);
                    
                    }else{
                    
                        complexPolyFusionCheck(complexMonomerContainer[i-1]);
                    
                    }
        
                    complexMonomerContainer[i-1].waitingTime = 0;
            
            }

        }else{


            while(complexCrowdantContainer[i-nComplexMonomers-1].waitingTime<dt){

                complexCrowdantDiffuse(complexCrowdantContainer[i-nComplexMonomers-1]);
                complexCrowdantContainer[i-nComplexMonomers-1].waitingTime+=-log(random)/crowdantDiffusionRate;
            }
                complexCrowdantContainer[i-nComplexMonomers-1].waitingTime-=dt;

        }
    }
}

// function for take a step
void System ::makeMove(int &x, int &y, int &newX, int &newY){

    int tempShift = lattice[x][y];
    lattice[x][y] = 0;
    lattice[newX][newY] = tempShift;
    x = newX;
    y = newY;
}

void System ::nextMove(int &oldX, int &oldY, int &diff, int &newX, int &newY, int& move){

    int tempRandom = random * 4;
    move = tempRandom;
    switch (tempRandom)
    {

    case 0:

        newX = (NX + (oldX + (diff)) % NX) % NX;
        newY = oldY;

        break;

    case 1:
        newX = (NX + (oldX - (diff)) % NX) % NX;
        newY = oldY;

        break;

    case 2:
        newY = (NY + (oldY + (diff)) % NY) % NY;
        newX = oldX;

        break;

    case 3:
        newY = (NY + (oldY - (diff)) % NY) % NY;
        newX = oldX;

        break;

    default:
        std ::cout << "Floor failed" << std ::endl;
    }
}


/*
void System ::polymerDiffuse()
{

    float normx = 0, normy = 0;
    int head1x = 0, head2x = 0, head1y = 0, head2y = 0, tempx = 0, tempy = 0, flag = 0;
    for (auto i : System ::polymerContainer)
    {

        float diffusion = monomerContainer[i.second.head1 - 1].diff / i.second.length;
        // std :: cout << diffusion << std :: endl;
        //  transverse
        if (diffusion <= 1)
        {

            if (random <= diffusion)
            {

                head1x = monomerContainer[i.second.head1 - 1].x;
                head1y = monomerContainer[i.second.head1 - 1].y;

                head2x = monomerContainer[i.second.head2 - 1].x;
                head2y = monomerContainer[i.second.head2 - 1].y;

                vectorComp(head1x, head1y, head2x, head2y, normx, normy);
                float distance = pow(pow(double(head1x - head2x), 2) + pow(double(head1y - head2y), 2), .5);

                if ((distance > (i.second.length - 1)) || (distance < (i.second.length - 1)))
                {
                    // broken=1;
                    if (normx != 0)
                    {

                        normx *= -1;
                    }
                    else if (normy != 0)
                    {

                        normy *= -1;
                    }
                }
                else
                {
                    // broken=0;
                }

                // std :: cout << normx << " " << normy << "yes";

                if (random < 0.5)
                {

                    tempx = head2x + normx;
                    tempy = head2y + normy;
                    flag = 1;
                }
                else
                {

                    normx *= -1;
                    normy *= -1;
                    tempx = head1x + normx;
                    tempy = head1y + normy;

                    flag = 2;
                }

                pbcX(tempx);
                pbcY(tempy);
                // std :: cout << head2x << " " << head2y << "\n";

                // std :: cout << tempx << " " << tempy << "yes";
                if (isVacant(tempx, tempy))
                {

                    for (auto j : i.second.seq)
                    {

                        // std :: cout << monomerContainer[j-1].x << " " << monomerContainer[j-1].y << std :: endl;
                        tempx = monomerContainer[j - 1].x + normx;
                        tempy = monomerContainer[j - 1].y + normy;

                        pbcX(tempx);
                        pbcY(tempy);

                        monomerContainer[j - 1].x = tempx;
                        monomerContainer[j - 1].y = tempy;

                        // std :: cout << monomerContainer[j-1].x << " " << monomerContainer[j-1].y << std :: endl;
                    }

                    if (flag == 1)
                    {

                        tempx = monomerContainer[i.second.head1 - 1].x - normx;
                        tempy = monomerContainer[i.second.head1 - 1].y - normy;

                        pbcX(tempx);
                        pbcY(tempy);

                        // std :: cout << tempx << " " << tempy << " " << lattice[tempx][tempy] << std :: endl;
                        lattice[tempx][tempy] = 0;
                        // std :: cout << tempx << " " << tempy << " " << lattice[tempx][tempy] << std :: endl;

                        tempx = monomerContainer[i.second.head2 - 1].x;
                        tempy = monomerContainer[i.second.head2 - 1].y;
                        lattice[tempx][tempy] = 1;
                    }
                    else if (flag == 2)
                    {

                        tempx = monomerContainer[i.second.head2 - 1].x - normx;
                        tempy = monomerContainer[i.second.head2 - 1].y - normy;

                        pbcX(tempx);
                        pbcY(tempy);

                        // std :: cout << tempx << " " << tempy << " " << lattice[tempx][tempy] << std :: endl;
                        lattice[tempx][tempy] = 0;
                        //  std :: cout << tempx << " " << tempy << " " << lattice[tempx][tempy] << std :: endl;

                        tempx = monomerContainer[i.second.head1 - 1].x;
                        tempy = monomerContainer[i.second.head1 - 1].y;
                        lattice[tempx][tempy] = 1;
                    }
                }
            }
            else
            {

                // No diffusion
            }
        }
        else
        {

            head1x = monomerContainer[i.second.head1 - 1].x;
            head1y = monomerContainer[i.second.head1 - 1].y;

            head2x = monomerContainer[i.second.head2 - 1].x;
            head2y = monomerContainer[i.second.head2 - 1].y;

            vectorComp(head1x, head1y, head2x, head2y, normx, normy);
            float distance = pow(pow(double(head1x - head2x), 2) + pow(double(head1y - head2y), 2), .5);

            if ((distance > (i.second.length - 1)) || (distance < (i.second.length - 1)))
            {

                if (normx != 0)
                {

                    normx *= -1;
                }
                else if (normy != 0)
                {

                    normy *= -1;
                }
            }

            if (random < 0.5)
            {

                tempx = head2x;
                tempy = head2y;
                flag = 1;
            }
            else
            {

                normx *= -1;
                normy *= -1;
                tempx = head1x;
                tempy = head1y;

                flag = 2;
            }

            int count = 0;
            for (auto k = 1; k <= diffusion; k++)
            {

                tempx += normx;
                tempy += normy;

                pbcX(tempx);
                pbcY(tempy);

                if (isVacant(tempx, tempy))
                {
                    count += 1;
                }
            }
            if (count == diffusion)
            {

                for (auto j : i.second.seq)
                {

                    // std :: cout << monomerContainer[j-1].x << " " << monomerContainer[j-1].y << std :: endl;

                    tempx = monomerContainer[j - 1].x + normx * diffusion;
                    tempy = monomerContainer[j - 1].y + normy * diffusion;

                    pbcX(tempx);
                    pbcY(tempy);

                    monomerContainer[j - 1].x = tempx;
                    monomerContainer[j - 1].y = tempy;

                    // std :: cout << monomerContainer[j-1].x << " " << monomerContainer[j-1].y << std :: endl;
                }

                head1x = monomerContainer[i.second.head1 - 1].x;
                head1y = monomerContainer[i.second.head1 - 1].y;

                head2x = monomerContainer[i.second.head2 - 1].x;
                head2y = monomerContainer[i.second.head2 - 1].y;

                if (flag == 1)
                {

                    auto it1 = i.second.seq.begin();
                    auto it2 = i.second.seq.rbegin();

                    for (int k = 1; k <= diffusion; k++)
                    {

                        tempx = monomerContainer[*it1 - 1].x - normx * k;
                        tempy = monomerContainer[*it1 - 1].y - normy * k;
                        pbcX(tempx);
                        pbcY(tempy);
                        lattice[tempx][tempy] = 0;

                        tempx = monomerContainer[*it2 - 1].x - normx * (k - 1);
                        tempy = monomerContainer[*it2 - 1].y - normy * (k - 1);
                        pbcX(tempx);
                        pbcY(tempy);
                        lattice[tempx][tempy] = 1;
                    }
                }
                else if (flag == 2)
                {

                    auto it1 = i.second.seq.rbegin();
                    auto it2 = i.second.seq.begin();

                    for (int k = 1; k <= diffusion; k++)
                    {

                        tempx = monomerContainer[*it1 - 1].x - normx * k;
                        tempy = monomerContainer[*it1 - 1].y - normy * k;
                        pbcX(tempx);
                        pbcY(tempy);
                        lattice[tempx][tempy] = 0;

                        tempx = monomerContainer[*it2 - 1].x - normx * (k - 1);
                        tempy = monomerContainer[*it2 - 1].y - normy * (k - 1);
                        pbcX(tempx);
                        pbcY(tempy);
                        lattice[tempx][tempy] = 1;
                    }
                }
            }
        }
    }
}
*/