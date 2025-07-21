// print.cpp
#include "class.hpp" 

// function for prinity details for monomers and crowdants
void System :: printEntity(std :: ofstream& monomerFile,std :: ofstream& crowdantFile ){

        int tempc = 0;
        monomerFile << "#start\t" << float(timeInternal) << std :: endl;
        for (complexMonomer& i: complexMonomerContainer){

            monomerFile  << std :: setw(10) << i.name << "\t" << std :: setw(5) << i.x << "\t" << std :: setw(5) << i.y << "\t" << std :: setw(3)<< "\t"<< i.mono << std :: setw(3) << "\t"  << std :: setw(3) << i.head << "\t" << std :: setw(3) << i.size << "\t" << std :: setw(6) << i.PID << std :: endl;                     

            tempc+=1;
            
            if(tempc>=1000){
                //std :: cout << "exit" << std :: endl;
                //break;
             }

        }
        monomerFile <<"#end"<< std :: endl;


        
        // for(crowdant& i : crowdantContainer){
        //    crowdantFile << std :: setw(10) << i.name << "\t" << std :: setw(5) << i.x << "\t" << std :: setw(5) << i.y << "\t" << std :: setw(3) << i.diff << "\t" << std :: setw(3) << 0  << std :: endl;
        
        //}
        

        crowdantFile << "#start\t" << float(timeInternal) << std :: endl;

        tempc = 0;
        for(complexCrowdant& i : complexCrowdantContainer){
            

            //std :: cout << i.name << i.size << std :: endl;
        crowdantFile << std :: setw(10) << i.name << "\t" << std :: setw(5) << (NX+(i.x)%NX)%NX << "\t" << std :: setw(5) << (NY+(i.y)%NY)%NY << "\t" << std :: setw(3)<< i.diff << "\t" << std :: setw(3) << i.size << std :: endl;


            // for(int m = -i.size;m<=i.size;++m){
            //     for(int n = -i.size;n<=i.size;++n){
            //                    crowdantFile << std :: setw(10) << i.name << "\t" << std :: setw(5) << (NX+(i.x+m)%NX)%NX << "\t" << std :: setw(5) << (NY+(i.y+n)%NY)%NY << "\t" << std :: setw(3)<< i.diff << "\t" << std :: setw(3) << i.size << std :: endl;
                                                                
            //    std :: cout << i.name << "\t" << i.x+m << "\t" << i.y+n << "\t" << 2 << "\t" << -1 << std :: endl;
             
            //     }                    
            // }
        
            tempc+=1;
            
            if(tempc>=1000){
                //std :: cout << "exit" << std :: endl;
              //  break;
             }
        }
        crowdantFile <<"#end"<< std :: endl;
        
}

    
// function for printing the lattice in the sparse form
void System :: printLattice(std :: ofstream& outfile){

        int temp=0;
        outfile << "#Time\t" << float(timeInternal) << std :: endl;

        for(int row=0;row<NX;row++){
            for(int col=0;col<NY;col++){
                if(lattice[row][col]>0){
                    outfile << row << "\t" << col << "\t" << lattice[row][col] << "\n";
                    temp+=1;
                }
            }
        }	
        
        std :: cout << monomerMass+crowdantMass+polymerMass << " " << temp << std :: endl;
        outfile << std :: endl;
        if(temp!=monomerMass+crowdantMass+polymerMass){
             std :: cout << nComplexMonomers << std :: endl;
             std :: cout << monomerMass << " " << crowdantMass << " " << polymerMass << std :: endl;
             std :: cout << temp << std :: endl;
             std :: cout << "exit 1" << std :: endl;
            exit(-1);}
}

void System :: printComplexPolymers(std :: ofstream& outfile){

        //tempPolymer=0;
        outfile << "#start\t" << float(timeInternal) << std :: endl;        


        for(auto& i : polymerContainer){                                

            //std :: cout << i.first << std :: endl;
            //std :: cout << i.second.length << std :: endl;

            
            float nx=0, ny=0;
            
            int h1x = complexMonomerContainer[i.second.head1-1].x;
            int h1y = complexMonomerContainer[i.second.head1-1].y;

            int h2x = complexMonomerContainer[i.second.head2-1].x;
            int h2y = complexMonomerContainer[i.second.head2-1].y;

          
            vectorComp(h1x,h1y,h2x,h2y,nx,ny);

            outfile << ">\t" <<i.second.PID << "\t" << i.second.length << "\t";
            //tempPolymer+=i.second.length;
            for(int& j : i.second.seq){
                outfile << j << "\t";
            }                                
                    
            outfile << nx << "\t" << ny << std :: endl;
        } 
        outfile << "#end" << std :: endl;        

}



// function for printing simulation parameters in a file
void System :: printLogFileHeader(std :: ofstream& outfile,std :: vector<std :: vector<float>> temp, int monomerSize, int totalTime, float iGridSpacing){

        outfile << "#Monomers               " << nComplexMonomers << "\t" << System :: monodiff  << "\t" << monomerSize << "\t" << monomerDiffusionRate << "\t" << monomerStepLength <<  std :: endl 
        << "#size X                 " << NX                   << std :: endl 
        << "#size Y                 " << NY                   << std :: endl
        << "#Grid spacing           " << iGridSpacing         << std :: endl 
        << "#Distribution Pattern   " << distributionPattern  << std :: endl
        << "#Total Time             " << totalTime*dt         << std :: endl
        << "#dt                     " << dt                   << std :: endl
        << "#reactionDt             " << reactionDt           << std :: endl
        << "#Monomer Binding Rate   " << monoBindRate         << std :: endl
        << "#Binding Rate           " << bindRate             << std :: endl
        << "#Unbinding Rate Nucleus " << unbindRateNucleus    << std :: endl
        << "#Unbinding Rate Polymer " << unbindRatePolymer    << std :: endl
        << "#Pre-seeds              " << preseededPolymers    << std :: endl
        << "#Pre-seed-size          " << preseededPolymerSize << std :: endl
        << "#Nucleus-size           " << nucleusSize          << std :: endl;
        for(uint i = 0; i < temp.size(); i++){  
        outfile << "#crowdant               " <<  i+1 << "\t" << temp[i][0] << "\t" << temp[i][1] << "\t" << temp[i][2] << "\t" << crowdantDiffusionRate << std :: endl;
        } 
        outfile << "#####################" 
        << std :: endl;          
        outfile << fmt << float(timeInternal) << "\t" << std :: setw(10) << 0                 << "\t" << std :: setw(5)  << 0               << "\t" << std :: setw(10) << polymerMass << "\t" << std :: setw(10) << monomerMass << "\t" << std :: setw(10) << nCrowdants << std :: endl;
           
}

// function for printing time in a file 
void System :: printLogs(std :: ofstream& outfile){
          
        outfile << fmt << float(timeInternal) << "\t" << std :: setw(10) << duration1.count() << "\t" << std :: setw(5) << duration2.count() << "\t" << std :: setw(10) << polymerMass << "\t" << std :: setw(10) << monomerMass << "\t" << std :: setw(10)  << nCrowdants << std :: endl;
    
}
