#include <fstream>
#include "../include/gyro.h"

namespace ouin
{

inline void wtf(const std::string& filename , const::std::string& in){
    std::fstream outfile;
    outfile.open(filename, std::fstream::in |std::fstream::out | std::fstream::app);
    if(!outfile){
        std::cout<<"creating new file..";
        outfile.open(filename, std::fstream::in |std::fstream::out | std::fstream::trunc) ;
        outfile <<"\n";
    }

    /*for(int z=0; z<gyro::icntl.size(); z++){
        outfile << gyro::icntl[z] << "\n" ;
    }*/
    outfile << in <<"\n" ;
    outfile.close();
}

}