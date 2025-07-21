#include "class.hpp"
    
// class constructor  
System :: System(int iNX, int iNY, float idt, int  iGridSpacing){
            
        NX              = iNX;                               // XMAX 
        NY              = iNY;                               // YMAX 
        dt              = idt;                               // integration time 
        gridSpacing     = iGridSpacing;
                        
        // setting lattice to be zero 
        for (int row=0;row<NX;row++){
        lattice.push_back(std :: vector<int>(NY,0));
            }            
}

void System :: checkIDS(){

    setIDS();

    for(int i =0;i<preseededPolymers;i++){
       // std :: cout << *IDS.begin() << std :: endl;
        IDS.erase(IDS.begin());
      //  std :: cout << *IDS.begin() << std :: endl;

    }

        startTime = std::chrono::high_resolution_clock::now();
        stepTime  = startTime;    
}

void System :: timeOn(){

        System :: stepTime = std::chrono::high_resolution_clock::now();

}

void System :: timeOff(){

        System :: currentTime = std::chrono::high_resolution_clock::now();

        duration1 = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime);
        duration2 = std::chrono::duration_cast<std::chrono::seconds>(currentTime - stepTime);


}
