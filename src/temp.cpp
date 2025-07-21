// temp

#include <iostream>
#include <cmath>
#include <vector> 
#include <random>
static std::random_device rd;
static std::mt19937 engine(rd());
static std::uniform_real_distribution<double> random_domain(0.0 , 1.0);
#define random random_domain(engine)



//#include <pair>
using namespace std;



void vectorComp(int& h1x, int& h1y, int& h2x, int& h2y, float& normx, float& normy){


        float dx2  = pow(h2x-h1x,2);
        float dy2  = pow(h2y-h1y,2);

        float mag = pow(dx2+dy2,.5);

        normx     = (h2x-h1x)/mag;
        normy     = (h2y-h1y)/mag;


}


void complexMonomerBindingOrientation(int ix, int iy, int jx, int jy, int  isize, int jsize, int &newx, int &newy){

    /*
    inputs to the function  
    
    ix is the x coordinate of monomer i 
    jx is the x coordinate of monomer j

    iy is the y coordinate of monomer i
    jy is the y coordinate of monomer j

    isize is the size of monomer i
    jsize is the size of monomer j  
    */

    /*
    outputs to the function 

    newx new x coordinate of monomer j 
    newy new y coordinate of monomer j
    */

    float compx,compy;
    vectorComp(ix,iy,jx,jy,compx,compy);
    
    std :: cout << "orientation  " << compx << " " << compy << " ";

    int xdist  = ix - jx;                 
    int ydist  = iy - jy;
    
    int absxdist = abs(xdist);
    int absydist = abs(ydist);

    cout << xdist << " " << ydist << endl;
    cout << absxdist << " " << absydist << endl;


    if(fabs(compx+compy)==1.0){
        newx = ix + (isize + jsize + 1)*compx;
        newy = iy + (isize + jsize + 1)*compy;

        std :: cout << 0 << " " << fabs(compx+compy) << std :: endl;
        return;
    }


    if(absxdist<absydist){
        if(ydist>=0){
        
            newx = ix; 
            newy = iy - isize - jsize - 1; 
            std :: cout << 1 << std :: endl;
        
        }else if(ydist<0){
        
            newx = ix;
            newy = iy + isize + jsize + 1;
            std :: cout << 2 << std :: endl;

        }
    }else if(absxdist>absydist){
        if(xdist>=0){
        
            newx = ix - isize - jsize - 1;
            newy = iy; 
            std :: cout << 3 << std :: endl;

        
        }else if(xdist<0){
        
            newx = ix + isize + jsize + 1;
            newy = iy; 
            std :: cout << 4 << std :: endl;

        }
    }else if(absxdist==absydist){

        if(random<.5){
            if(ydist>=0){
            
                newx = ix; 
                newy = iy - isize - jsize - 1; 
                std :: cout << 5 << std :: endl;

            
            }else if(ydist<0){
            
                newx = ix;
                newy = iy + isize + jsize + 1;
                std :: cout << 6 << std :: endl;

            }
        }else{
            if(xdist>=0){
            
                newx = ix - isize - jsize - 1;
                newy = iy; 
                std :: cout << 7 << std :: endl;

            }else if(xdist<0){

                newx = ix + isize + jsize + 1;
                newy = iy; 
                std :: cout << 8 << std :: endl;

            }
        }
    }
    return;
}



int main(){

    int x=3,y=3,size=1,tx,ty;
    complexMonomerBindingOrientation(4,15,6,17,1,1,tx,ty);
    std :: cout << tx <<  " " << ty << std :: endl;

}