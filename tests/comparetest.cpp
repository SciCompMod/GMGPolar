#include <gtest/gtest.h>
#include <fstream>
#include <sched.h>
#include <time.h>
#include "../include/cmdline.h"
#include "../include/gmgpolar.h"

class SimpleTest: public ::testing::Test{
    protected:
    void SetUp() override{
        std::ifstream resource(Param::filename);
        std::string line;
        if(resource.is_open()){
            while(getline(resource,line)){
                input.push_back(line);
            }
        }
        initparam = 6;
        int initarr[initparam]= {Param::prob, Param::alpha_coeff, Param::beta_coeff, Param::nr_exp, Param::ntheta_exp, Param::mod_pk};
        data.resize(initparam);
        std::copy(initarr,initarr+initparam , data.begin());

        //Initializing default parameters
        gyro::init_params();
        gyro::icntl[Param::verbose]= 0;
        gyro::icntl[Param::debug]=0;
        gyro::icntl[Param::extrapolation]=0;
        gyro::icntl[Param::DirBC_Interior]=1;
        gyro::icntl[Param::check_error]=1;
        gyro::dcntl[Param::R0]= 1e-5;
        gyro::f_grid_r = "";
        gyro::f_grid_theta= "";
        gyro::f_sol_in= "";
        gyro::f_sol_out="";
        gyro::icntl[Param::nr_exp] = 4;
        gyro::icntl[Param::ntheta_exp]= 4 ;
        gyro::icntl[Param::fac_ani]=3 ;
        gyro::select_functions_class(gyro::icntl[Param::alpha_coeff], gyro::icntl[Param::beta_coeff], gyro::icntl[Param::mod_pk], gyro::icntl[Param::prob]);

    }

    int initparam;
    std::vector<std::string> input;
    std::vector<int> data;
};

TEST_F(SimpleTest, Initialize){
    for(int z=0 ;z<initparam; z++){
        ASSERT_EQ(gyro::icntl[data[z]],std::stoi(input[z])) <<"No use testing if initial data is different";
    }
}

TEST_F(SimpleTest, finestgrid){
    gmgpolar gmgtest;
    gmgtest.create_grid_polar();//only the finest grid is now created
    EXPECT_EQ(gmgtest.v_level.size(),1);

    int finest_nodes = gmgtest.v_level[0]->nr * gmgtest.v_level[0]->ntheta;
    
    EXPECT_EQ( finest_nodes , std::stoi( input[initparam+1] ) ) ;
    //EXPECT_EQ( gmgtest.v_level[0]->nr, std::stoi( input[initparam+2])); Segmentation fault ???
}

TEST_F(SimpleTest, Multigrid){
    gmgpolar gmgtest2;
    gmgtest2.create_grid_polar();
    gmgtest2.polar_multigrid();

    EXPECT_EQ(gmgtest2.levels, std::stoi( input[initparam+4 ] ) );
    //TODO: iteration number
}