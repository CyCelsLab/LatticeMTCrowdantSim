#pragma once
/*
Aman Soni, 
Indian Institute of Science Education and Research, Pune 
c++ script for simulating the effect of macromolecular on polymerisation kinetics 
*/

#include<iostream>
#include<cstdio>
#include<chrono>
#include<fstream>
#include<string>
#include<list>
#include<set>

#include<unordered_map>
#include<iomanip>
#include<cmath>
#include<vector>
#include<algorithm>

// formatting snippet 
#define fmt  std :: setw(14) << std :: setprecision(5) << std :: scientific


class System{

    public:
      
        /* structure of monomer, crowdants and polymer*/
      

        // structure of polymer
        struct polymer{
        
            int  PID       = 0;             //  name of the polmer
            int  head1     = 0;             //  head1 of the polymer 
            int  head2     = 0;             //  head2 of the polymer   
            int  length    = 0;             //  length of the polymer    
            
            std :: list<int> seq = {};      //  polymer sequence   
            
            int  h1x       = -1;
            int  h1y       = -1;
            int  h2x       = -1;
            int  h2y       = -1; 

        };        
//-------------------------------------------------------------------------------------------------            

        // structur of complex crowdant
        struct complexCrowdant{

            int name   =  0;
            int x;
            int y;
            int size =  0;
            int diff   =  0; 
            float waitingTime = 0.0;
            int mass   = 0;
           // int diffTime = 1;

        };

//-------------------------------------------------------------------------------------------------

        // structure for monomers 
        struct complexMonomer{
    
            int  name =  0;           //   name of the monomer
            int  x    = -1;           //   x position  
            int  y    = -1;           //   y position
            int  mono =  true;        //   true if not in any polymer  
            int  head =  false;       //   true if head of ant polymer  
            int  PID  = -1;           //   not -1 if part of any polymer  
            int  size =  0;
            float waitingTime = 0.0;
            float reactionDt  = 0.0;
           //int diffTime = 1;

        };
//-------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------            

        // function for imposing PBC in x direction
        void pbcX(int &x){
            x = (NX+(x)%NX)%NX;
        }

        // function for imposing PBC in y direction
        void pbcY(int &y){
            y = (NY+(y)%NY)%NY;
        }


        // class constructor 
        
        //---------Constructor.cpp----------// 
        System(int iNX, int iNY, float idt, int iGridSpacing);

        // function to start timer 
        void timeOn();
        // function of end timer 
        void timeOff();
        
        //---------public methods-----//
        
        
        //-------populating.cpp------//
        
        bool checkBlockInLattice(int& x, int& y, int& size);
  
        // function for filling the lattice with a specific value
        void fillBlockLattice(int& x, int& y, int& size, int value);

        // function to check if block of size x size is empty
        bool checkEmptyBlock(int& x, int& y, int& size);

        // function to check if block can be moved or not
        bool checkCanMoveBlock(int& x, int& y, int& size, int name);

        // taking complex crowdant in the system 
        void putComplexCrowdants(int icomplexCrowdantName, int size, int stepLength);

        // taking complex monomer in the system
        void putComplexMonomers(int iComplexMonomerName, int size, int stepLength);

        // generating polymer IDS 
        void checkIDS();


        //----print.cpp----//
        
        // function for prinity details for monomers and crowdants
        void printEntity(std :: ofstream& monomerFile,std :: ofstream& crowdantFile);
        
        // function for printing the lattice in the sparse form
        void printLattice(std :: ofstream& outfile);
        
        // function of printing polymer in a file
        void printComplexPolymers(std :: ofstream& outfile);
        
        // function for printing simulation parameters in a file
        void printLogFileHeader(std :: ofstream& outfile,std :: vector<std :: vector<float>> temp, int monomerSize, int totalTime, float iGridSpacing);
        
        // function for printing time in a file 
        void printLogs(std :: ofstream& outfile);
        

        
        //-----diffuse.cpp-----//

        // function for evolving system in diffusion and reaction 
        void diffuseandReact();
        
        //-----polymerise.cpp-----//
                
        // function for polymersing complex monomers
        void complexPolymerise();


        // function for providing the unit direction of the vector.
        void vectorComp(int& h1x, int& h1y, int& h2x, int& h2y, float& normx, float& normy){


                float dx2  = pow(h2x-h1x,2);
                float dy2  = pow(h2y-h1y,2);

                float mag = pow(dx2+dy2,.5);

                normx     = (h2x-h1x)/mag;
                normy     = (h2y-h1y)/mag;


        }

        void checkDepolyType(int condition){
            if (condition==1){
                depolType=1;
            }else{
                depolType=2;
            }

        }

        int getDpolType(){
            if(depolType==1){return (1);}
            else{return(2);}
        }



        //-----depolymerise.cpp-----//
        
        // function for depolymerising polymers 
        void depolymerise();

        // function for depolymerising at only one end
        void singleDepolymerise();
        //-----test.cpp-----//
        // testing code in testing and debugging 
        void test();




        //-----preseed.cpp-----//

        // for preseeding the system
        void preseeds(int l, int length, int size);


        // function for maintaining internal clock  
        void shiftTime(){
        
                currentTime = std::chrono::high_resolution_clock::now();
                timeInternal+=(dt);
                iterations+=1;
        }

        // setting the possible polymer ids 
        void setIDS(int k){

            for(int i=0;i<(k);i++){   
                IDS.insert(i+1);
            }

        }


        //----Kinetics.cpp-----//
        // function for adding kinestics in the model

        void putKinetics(float imonoBindRate, float ibindRate, float iunbindRateNucleus , float iunbindRatePolymer, float reactionDt,int inucleusSize);

        void populateLattice(int iTotalComplexCrowdants, int csize, int cdiffStepLength, int iTotalComplexMonomers, int msize, int mdiffStepLength, int ipreseeds, int ipreseedSize, int iDistributionPattern);


        std :: vector<long long> monoHopCounter   = {0,0,0,0};
        std :: vector<long long> monoTotalCounter = {0,0,0,0};
        std :: vector<long long> crowTotalCounter = {0,0,0,0};
        std :: vector<long long> crowHopCounter   = {0,0,0,0};

        void equillibrate(int iter);    


    private:
    
        //Private variables//
                
        // lattice (NX,NY) matrix
        std :: vector<std :: vector<int>> lattice;
        
        // hash table for storing polymer details
        std :: unordered_map<int,polymer> polymerContainer;
                    

        // container for complex crowdants
        std :: vector<complexCrowdant> complexCrowdantContainer;

        // vector contains the details of crowdant jump intervals 
        std :: vector<std :: vector<int>> crowdantJumpInterval;

        // container for complex monomers 
        std :: vector<complexMonomer> complexMonomerContainer;
        //std :: vector<std :: vector<int>> monomerJumpInterval;
    
        // containeTimer for storing available polymer ids
        std :: set<int> IDS;       



        // total number of monomers 
        int totalMonomers   = 0;
        
        // size of the matrix
        int NX       = 0;
        int NY       = 0;
        int gridSpacing = 1;
                
        // crowdant mass

        // number of crowdants present in the system 
        int nCrowdants = 0;

        // number of complex crowdants present in the system
        int nComplexCrowdants = 0;

        // number of complex monomers pressent in the system 
        int nComplexMonomers = 0; 

        // integration time 
        float           dt=0;
        float           reactionDt=0;
        
        // monomer binding rate and probability 
        float monoBindRate;
        float monoBindProb;
        
        // binding rate and proability 
        float bindRate;
        float bindProb;
        
        // unbinding rate and probability for polymer 
        float unbindRatePolymer;
        float unbindProbPolymer;

        // unbinding rate and probability for nucleus 

        float unbindRateNucleus;
        float unbindProbNucleus;
        
        // time step counter 
        float timeInternal = 0;
        int   iterations   = 0;
        
        // number of polymer and monomers 
        int  monomerMass  = 0;
        int  polymerMass  = 0;
        int  crowdantMass = 0;
        int  totalEntity  = 0;

        int preseededPolymers=0;
        int preseededPolymerSize=0;
        int nucleusSize=2;

        //int monomerStepLength  = 1;
        //int crowdantStepLength = 1;

        float monomerDiffusionRate  = 0;
        float crowdantDiffusionRate = 0; 

        float monodiff = 0;
        float crowdiff = 0;

        float depolType=2;


        int monomerStepLength  = 0;
        int crowdantStepLength = 0;

        int distributionPattern = 0;
        int equllibriationTime  = 0;

        // number of tatal enetities present in the system (monomers+polymers+crowdants)
        int totalMass   = 0;

        // neighbour vector
            std :: vector<std :: pair<int,int>> neighbourVector = {{0,1},{0,-1},{-1,0},{1,0}};


        // Monomer index for random iterations 

        std :: vector<int> randomEntity;
        
        // high resolution clock for timming the time lapsed          
        std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now(); 
        std::chrono::time_point<std::chrono::high_resolution_clock> stepTime  = std::chrono::high_resolution_clock::now();
        std::chrono::time_point<std::chrono::high_resolution_clock> currentTime;
              
        std::chrono:: seconds  duration1;      
        std::chrono:: seconds  duration2;      





        //-----populating.cpp-----//
        // function for adding crowdants
        
        // setting the possible polymer ids 
        void setIDS();
            
        //-----polymerise.cpp-----//
        // mono-mono fusion     

        void complexMonoPolyBindingOrientation(complexMonomer& m , complexMonomer& n, int& tempx,int& tempy, int& head);

        // searching neighbours for complex monomers.
        void neighbours(int x, int y, int size, std :: set<int> &tempset);

        void calculateMonoAxis(complexMonomer &i, int newx, int newy, polymer &p);

        void complexMonoFusion(complexMonomer &i, complexMonomer &j, int newx, int newy);

        void complexMonoFusion_(complexMonomer &i, complexMonomer &j, int newx, int newy);

        void complexMonoMonoBindingOrientation(int ix, int iy, int jx, int jy, int  isize, int jsize, int &newx, int &newy);

        void complexMonoMonoBindingOrientation_(int ix, int iy, int jx, int jy, int isize, int jsize, int &newx, int &newy);


        void complexMonoPolyBindingOrientation_(complexMonomer &m, complexMonomer &n, int &tempx, int &tempy, int &head);

        void complexMonoPolyFusion(complexMonomer &n, complexMonomer &m, int head, int PID, int tempx, int tempy);

        void complexMonoFusionCheck(complexMonomer& i);

        void complexPolyFusionCheck(complexMonomer& i);


        // to know which end of the polymer.
        void whichHead(int& head1, int& head2, int& name, int& head);
   
        // functhead1ion to check if lattice site is empty
        bool isVacant(int& i, int& j){
            bool vacant = lattice[i][j]==0;    
            return vacant;
        }
        
        //-----Diffuse.cpp-----//
        
        // function for diffusing crowdants 
        void crowdantDiffuse();
        
        // function for diffusing monomers 
        void monomerDiffuse();

        // function of diffusing complex crowdants
        void complexCrowdantDiffuse(complexCrowdant &i);

        // function of diffusing complex monomer
        void complexMonomerDiffuse(complexMonomer &i);

        // function of diffusing polymer 
        void polymerDiffuse();


        // function for getting next position for diffusible entity
        void nextMove(int& oldX, int& oldY, int& diff, int& newX, int& newY, int& move);

        // function for take a step 
        void makeMove(int& x, int& y, int& newX, int& newY);

};
