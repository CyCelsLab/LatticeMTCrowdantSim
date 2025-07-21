#include "class.hpp"
#include "random.hpp"

// function for depolymerising polymers 
void System :: depolymerise(){


    for(auto iter = polymerContainer.begin(); iter!=polymerContainer.end();){

            
        auto it = iter++;


        if(it->second.length==2){
            if(random<=unbindProbNucleus){
                
                int PPID = it->first;
                
                //int H1PID = monomerContainer[it->second.head1-1].PID;
                //int H2PID = monomerContainer[it->second.head2-1].PID;
                
                
                complexMonomerContainer[it->second.head2-1].mono =  1;
                complexMonomerContainer[it->second.head2-1].head =  1;
                complexMonomerContainer[it->second.head2-1].PID  = -1;
                complexMonomerContainer[it->second.head2-1].waitingTime = -log(random)/monomerDiffusionRate;
                complexMonomerContainer[it->second.head2-1].reactionDt  =  complexMonomerContainer[it->second.head2-1].waitingTime;


                complexMonomerContainer[it->second.head1-1].mono =  1;
                complexMonomerContainer[it->second.head1-1].head =  1;
                complexMonomerContainer[it->second.head1-1].PID  = -1;
                complexMonomerContainer[it->second.head1-1].waitingTime = -log(random)/monomerDiffusionRate;
                complexMonomerContainer[it->second.head1-1].reactionDt  =  complexMonomerContainer[it->second.head1-1].waitingTime;



                IDS.insert(PPID);
                
                
                polymerContainer.erase(it);
                
                int tempMass = pow(2*complexMonomerContainer[it->second.head2-1].size+1,2);
            
                polymerMass-=tempMass*2;
                monomerMass+=tempMass*2;                            
            }




        }else if(it->second.length>2){

            int temprandomint = random*2;
            
            float tempRate  = 0.00;
            if(it->second.length<=nucleusSize){
                tempRate = unbindProbNucleus;
            }else{
                tempRate = unbindProbPolymer;
            }


            if(temprandomint==0){
            
                //for head1
                if(random<tempRate){
                
                    complexMonomerContainer[it->second.head1-1].mono = 1;
                    complexMonomerContainer[it->second.head1-1].head = 0;
                    complexMonomerContainer[it->second.head1-1].PID = -1;
                    complexMonomerContainer[it->second.head1-1].waitingTime = -log(random)/monomerDiffusionRate;
                    complexMonomerContainer[it->second.head1-1].reactionDt  =  complexMonomerContainer[it->second.head1-1].waitingTime;
                    
                    it->second.seq.pop_front();
                    it->second.length-=1;
                    it->second.head1 = it->second.seq.front();

                    complexMonomerContainer[it->second.head1-1].head = 1;



                    int tempMass = pow(2*complexMonomerContainer[it->second.head2-1].size+1,2);            
                    polymerMass-=tempMass;
                    monomerMass+=tempMass;



                }
            }else{
                //for head2
                if(random<tempRate){
                
                    complexMonomerContainer[it->second.head2-1].mono = 1;
                    complexMonomerContainer[it->second.head2-1].head = 0;
                    complexMonomerContainer[it->second.head2-1].PID = -1;
                    complexMonomerContainer[it->second.head2-1].waitingTime = -log(random)/monomerDiffusionRate;
                    complexMonomerContainer[it->second.head2-1].reactionDt  =  complexMonomerContainer[it->second.head2-1].waitingTime;


                    it->second.seq.pop_back();
                    it->second.length-=1;
                    it->second.head2 = it->second.seq.back();

                    complexMonomerContainer[it->second.head2-1].head = 1;

                    int tempMass = pow(2*complexMonomerContainer[it->second.head2-1].size+1,2);
            
                    polymerMass-=tempMass;
                    monomerMass+=tempMass;
                    
                }
            }    
        }

    }

}


void System :: singleDepolymerise(){

    for(auto iter = polymerContainer.begin(); iter!=polymerContainer.end();){

            
        auto it = iter++;

        if(it->second.length==2){

            if(random<=unbindProbNucleus){
                
                int PPID = it->first;
                
                //int H1PID = monomerContainer[it->second.head1-1].PID;
                //int H2PID = monomerContainer[it->second.head2-1].PID;
                
                
                complexMonomerContainer[it->second.head2-1].mono =  1;
                complexMonomerContainer[it->second.head2-1].head =  1;
                complexMonomerContainer[it->second.head2-1].PID  = -1;
                complexMonomerContainer[it->second.head2-1].waitingTime = -log(random)/monomerDiffusionRate;
                complexMonomerContainer[it->second.head2-1].reactionDt  =  complexMonomerContainer[it->second.head2-1].waitingTime;


                complexMonomerContainer[it->second.head1-1].mono =  1;
                complexMonomerContainer[it->second.head1-1].head =  1;
                complexMonomerContainer[it->second.head1-1].PID  = -1;
                complexMonomerContainer[it->second.head1-1].waitingTime = -log(random)/monomerDiffusionRate;
                complexMonomerContainer[it->second.head1-1].reactionDt  =  complexMonomerContainer[it->second.head1-1].waitingTime;



                IDS.insert(PPID);
                
                
                polymerContainer.erase(it);
                
                int tempMass = pow(2*complexMonomerContainer[it->second.head2-1].size+1,2);
            
                polymerMass-=tempMass*2;
                monomerMass+=tempMass*2;                            
            }




        }else if(it->second.length>2){

            float tempRate  = 0.00;
            if(it->second.length<=nucleusSize){
                tempRate = unbindProbNucleus;
            }else{
                tempRate = unbindProbPolymer;
            }

            //for head1
            if(random<tempRate){
            
                complexMonomerContainer[it->second.head1-1].mono = 1;
                complexMonomerContainer[it->second.head1-1].head = 0;
                complexMonomerContainer[it->second.head1-1].PID = -1;
                complexMonomerContainer[it->second.head1-1].waitingTime = -log(random)/monomerDiffusionRate;
                complexMonomerContainer[it->second.head1-1].reactionDt  =  complexMonomerContainer[it->second.head1-1].waitingTime;
                
                it->second.seq.pop_front();
                it->second.length-=1;
                it->second.head1 = it->second.seq.front();

                complexMonomerContainer[it->second.head1-1].head = 1;

                int tempMass = pow(2*complexMonomerContainer[it->second.head2-1].size+1,2);            
                polymerMass-=tempMass;
                monomerMass+=tempMass;
            }
        }
    }
}

