//test.cpp

#include "class.hpp"
#include "random.hpp"

// testing code in testing and debugging 
void System :: test(){

        polymer test;
        /*
        test.PID = 3;
        test.head1 = 7;
        test.head2 = 8;
        test.length = 2;
        test.seq = {7,8};
        polymerContainer.insert({3,test});

        test.PID = 2;
        test.head1 = 1;
        test.head2 = 2;
        test.length = 2;
        test.seq = {1,2};
        
        
        polymerContainer.insert({2,test});

        test.PID = 1;
        test.head1 = 4;
        test.head2 = 6;
        test.length = 3;
        test.seq = {4,5,6};

        polymerContainer.insert({1,test});

        //for(int i =0;i<3;i++) IDS.erase(IDS.begin());
        //for(auto i : IDS) cout << i << endl;

        monomer temp;

        for(int i=0;i<3;i++){
            
            temp.name = (i+1);
            temp.x    = 1;
            temp.y    = (i+1);                                     
            temp.diff = 3;
            temp.mono = 0;
            temp.PID  = 2;
            temp.head = 1;
            monomerContainer.push_back({temp});
            lattice[1][i+1] = 1;
        }
        
        monomerContainer[2].mono = 1;
        monomerContainer[2].head = 1;
        monomerContainer[2].PID  = -1;
                
        for(int i=3;i<6;i++){
            
            temp.name = (i+1);
            temp.x    = 2;
            temp.y    = (i+1);                                     
            temp.diff = 3;
            temp.mono = 0;size
            temp.head = 1;
            temp.PID  = 1;
            monomerContainer.p    for(auto i =0;i<iTotalPreSeeds;i++){

        

    }
ush_back(temp);
            lattice[2][i+1] = 1;
            
        }

        monomerContainer[4].head = 0;

        for(int i=6;i<8;i++){
            
            temp.name = (i+1);
            temp.x    = 3;
            temp.y    = (i+2);                                     
            temp.diff = 3;
            temp.mono = 0;
            temp.head = 1;
            temp.PID  = 3;
            monomerContainer.push_back(temp);
            lattice[5][i+1] = 1;
        }

        */

        test.PID = 1;
        test.head1 = 1;
        test.head2 = 4;
        test.length = 4;
        test.seq = {1,2,3,4};
        polymerContainer.insert({1,test});

        complexMonomer temp;

        for(int i=0;i<4;i++){
            
            temp.name = (i+1);
            temp.x    = (10);
            temp.y    = (10-i);

            temp.mono = 0;
            temp.PID  = 2;
            temp.head = 0;
            temp.size = 0; 
            complexMonomerContainer.push_back({temp});
            lattice[temp.x][temp.y] = 1;
        }

        complexMonomerContainer[1].head=1;
        complexMonomerContainer[2].head=1;
    

}

