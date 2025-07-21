#include "class.hpp"
#include<string>
#include<fstream> 
#include<filesystem>
#include <sstream>
#include <yaml-cpp/yaml.h>
#include <fstream>



int main(int argc, char** argv){


    // initialising required files 
    std :: ofstream entityFile,polymersFile,latticeFile, logFile, monomerFile, crowdantFile, encounterProb;    



    // opening required files
    std :: string src = "";
    if(argc<2){
        src.append("Data/1");
        if(std :: filesystem :: is_directory(src)==false){
            std :: filesystem :: create_directories(src);
        }
    }else{
        src.append( "Data/" + std :: string(argv[1]));
        if(std :: filesystem :: is_directory(src)==false){
            std :: filesystem :: create_directories(src);
        }
    }


    // Open the YAML file for parsing
    std::ifstream fin("config.yaml");
    if (!fin) {
        std::cerr << "Failed to open the YAML file.\n";
        return 1;
    }



//    // read the file config.txt 
//    std :: string input;
//    std :: ifstream MyReadFile("./config.txt");

    float   totalComplexMonomer, complexMonomerSize, complexMonomerStepLength;

    int     sizeX;
    int     sizeY;
    int     iGridSpacing=1;

    int     time;

    int     DistributionPattern = 0;
    float   dt          =.01;
    float   reactionDt  = dt;

    float   BRNucleus;
    float   BRPolymer;  
    float   UBRNucleus;
    float   UBRPolymer;

    int     NucleusSize;

    int     preseeds    = 0;
    int     preSeedSize = 0;

    int     crowdants   = 0; 

    int     printInterval  = 10;
    int     printLatttice  =  0;
    int     printPolymer   = 0;
    int     printMonoCrow  = 0;
    int     printLog       = 1;  
    int     depolType      = 2;
 
    std :: vector<std :: vector<float>> complexCrowdantsReaderList;
    std :: vector<std :: vector<float>> complexMonomersReaderList;



    YAML::Node config = YAML::Load(fin);


    // units  number 
    totalComplexMonomer       =  config["number-of-monomers"].as<float>();
        
    // [0,1,2,3,4] 
    complexMonomerSize        =  config["monomer-size"].as<float>();
    
    // 1e-7 meters
    complexMonomerStepLength  =  config["monomer-step-length"].as<float>();

    // always in lattice units (0,1,2,3..N)
    sizeX                     = config["sizeX"].as<int>();  
    
    // always in lattice units (0,1,2,3..N)    
    sizeY                     = config["sizeY"].as<int>();  

    // 1e-7 meters
    iGridSpacing        = config["grid-spacing"].as<float>();
    
    // seconds 
    time                = config["iterations"].as<float>();
    
    // [0,1,2] = [random, monomer first, crowdant first]
    DistributionPattern = config["distribution-pattern"].as<float>();
    
    // seconds 
    dt                  = config["integration-time"].as<float>();
    
    // seconds 
    reactionDt          = config["reaction-dt"].as<float>();    

    // 1/seconds		
    BRNucleus             = config["nucleus-binding-rate"].as<float>();
    // 1/seconds
    BRPolymer              = config["binding-rate"].as<float>();
    // 1/seconds
    UBRNucleus      = config["nucleus-unbinding-rate"].as<float>();
    // 1/seconds
    UBRPolymer      = config["unbinding-rate"].as<float>();

    // number 
    NucleusSize     = config["nucleus-size"].as<float>();

    // number 
    preseeds        = config["number-of-preseeds"].as<float>();
    
    // [0,1,2,3,4]
    preSeedSize     = config["presees-size"].as<float>();
    
    crowdants       = config["number-of-crowdant-type"].as<float>();
    
    // depolymerisation type 
    
    // 1 for single end depolymerisation 
    // 2 for depolymerisation at both ends 
    depolType       = config["depol-type"].as<int>();



    if(crowdants == 0){
        complexCrowdantsReaderList.push_back({0,0,0});

    }else{

        complexCrowdantsReaderList.push_back({});
	// number 
        complexCrowdantsReaderList[0].push_back(config["number-of-crowdants"].as<float>());
        	
	// [0,1,2,3,4] 
        complexCrowdantsReaderList[0].push_back(config["crowdant-size"].as<float>());

	// lattice units 
        complexCrowdantsReaderList[0].push_back(config["crowdant-step-length"].as<float>());    
    }


    printInterval    = config["print-interval"].as<int>();


    printLatttice = config["print-lattice"].as<int>();
    printPolymer  = config["print-polymer"].as<int>();
    printMonoCrow  = config["print-monomer-crowdant"].as<int>();
    printLog      = config["print-log"].as<int>();
 

    // opening files 
    latticeFile.open(src+"/Lattice.dat");
    polymersFile.open(src+"/Polymer.dat");
    logFile.open(src+"/Log.dat"); 
    monomerFile.open(src+"/Monomer.dat");
    crowdantFile.open(src+"/Crowdant.dat");
    encounterProb.open(src+"/Encounter.dat");

        
    // initialising the system populateLattice
    System tubulins(sizeX,sizeY,dt,iGridSpacing);

    // put kinetics parameters in the systems 
    tubulins.putKinetics(BRNucleus,BRPolymer,UBRNucleus,UBRPolymer,reactionDt,NucleusSize);


    // adding preseeds in the system
    tubulins.preseeds(preseeds,preSeedSize, complexMonomerSize);

    // 
    tubulins.populateLattice(complexCrowdantsReaderList[0][0],complexCrowdantsReaderList[0][1],complexCrowdantsReaderList[0][2],totalComplexMonomer,complexMonomerSize,complexMonomerStepLength,preseeds,preSeedSize, DistributionPattern);

    tubulins.checkDepolyType(depolType);

    // adding crowdants in the system 
//    for( auto i : complexCrowdantsReaderList){
//        std :: cout << i[0] << " " << i[1] << " " << i[2] <<  " " << i[3] << std :: endl;
//       // tubulins.putComplexCrowdants(i[0],i[1],i[2],complexCrowdantStepLength);
//    }

    
    // adding monomers in the system 
    //tubulins.putComplexMonomers(totalComplexMonomer,complexMonomerDiff,complexMonomerSize,complexMonomerStepLength,);  



    tubulins.equillibrate(500);



    // print entities in file 
    tubulins.printEntity(monomerFile,crowdantFile);

    // print occupied sites in lattice  
    // tubulins.printLattice(latticeFile);

     //print log file header 
    tubulins.printLogFileHeader(logFile,complexCrowdantsReaderList,complexMonomerSize,time,iGridSpacing);
    
    // polymer details if any 
    tubulins.printComplexPolymers(polymersFile);


    tubulins.checkIDS();



    // loop over number of iteration 
    for(int t = 1;t<=time;t++){


       if(t%10==0){
            std :: cout << "time " <<  t << std :: endl;
        }

        tubulins.shiftTime();

        tubulins.timeOn();

        tubulins.diffuseandReact();


        if(tubulins.getDpolType()==1){
            tubulins.singleDepolymerise();    
        }else{
            tubulins.depolymerise();    
        }
        tubulins.timeOff();

       
       // bookeeping 
       // --------------- // 
 
        if(t%(int(printInterval))==0){

            // print monomer and crowdants 
    	    if (printMonoCrow) tubulins.printEntity(monomerFile,crowdantFile);
            
            // print polymers 
            if (printPolymer) tubulins.printComplexPolymers(polymersFile);

            // print log file 
            if (printLog)     tubulins.printLogs(logFile);   

            // print lattice 
            if (printLatttice) tubulins.printLattice(latticeFile);
        }

    }

    tubulins.printComplexPolymers(polymersFile);
    tubulins.printLattice(latticeFile);
    
    for(auto i = 0;i < 4; i++){
        
        encounterProb << tubulins.monoHopCounter[i] << " " << tubulins.monoTotalCounter[i] << " " << tubulins.crowHopCounter[i] << " " << tubulins.crowTotalCounter[i] << std :: endl; 

    }

    // closing all the files 
    entityFile.close();
    latticeFile.close();
    polymersFile.close();
    logFile.close();    
// end of the file 
}
