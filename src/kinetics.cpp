// kinetics.cpp

#include "class.hpp"

void System :: putKinetics(float imonoBindRate, float ibindRate, float iunbindRateNucleus , float iunbindRatePolymer, float reactionDt, int inucleusSize){
    
        System :: reactionDt = reactionDt; 
        
        monoBindRate    = imonoBindRate;                     // monomer-monomer binding rate
        monoBindProb    = monoBindRate * reactionDt;         // monomer-monomer binding probability
        //monoBindProb    = 1;
                
        bindRate        = ibindRate;                         // polymer-monomer binding rate   
        bindProb        = bindRate * reactionDt;             // polymer-monomer binding prob   
        //bindProb        = 1;
        
        unbindRatePolymer      = iunbindRatePolymer;                       // unbinding rate    
        unbindProbPolymer      = unbindRatePolymer*reactionDt;             // unbinding probability    

        unbindRateNucleus      = iunbindRateNucleus;                       // unbinding rate    
        unbindProbNucleus      = unbindRateNucleus*reactionDt;             // unbinding probability    

        nucleusSize            = inucleusSize;

}
