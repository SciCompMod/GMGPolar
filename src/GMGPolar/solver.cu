#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

void GMGPolar::solve() {

//     const Level& level1 = levels_[0];
//     const PolarGrid& grid1 = level1.grid();
//     const PolarGrid& sgrid1 = level1.grid();
//     Vector<double> v1(grid1.numberOfNodes());
//     GPU_Vector<double> sv1(sgrid1.numberOfNodes());

//     const Level& level2 = levels_[1];
//     const PolarGrid& grid2 = level2.grid();
//     const PolarGrid& sgrid2 = level2.grid();
//     Vector<double> v2(grid2.numberOfNodes());
//     GPU_Vector<double> sv2(sgrid2.numberOfNodes());

//     for (int i = 0; i < v1.size(); i++){
//         v1[i] = i;
//     }
//     copyHostToDevice(v1,sv1);
//     for (int i = 0; i < v2.size(); i++){
//         v2[i] = i;
//     }
//     copyHostToDevice(v2,sv2);


//     GPU_Vector<double> sx1(sgrid2.numberOfNodes());
//     GPU_Vector<double> sy1(sgrid2.numberOfNodes());
//     GPU_Vector<double> sz1(sgrid2.numberOfNodes());
//     Vector<double> x1(sgrid2.numberOfNodes());
//     Vector<double> y1(sgrid2.numberOfNodes());
//     Vector<double> z1(sgrid2.numberOfNodes());
//     for (int i = 0; i < x1.size(); i++){
//         x1[i] = i;
//         y1[i] = 3*i;
//         z1[i] = 2*i;
//     }
//     copyHostToDevice(x1,sx1);
//     copyHostToDevice(y1,sy1);
//     copyHostToDevice(z1,sz1);



// auto start = std::chrono::high_resolution_clock::now();
//    level2.computeResidual(sx1, sy1, sz1);

//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Residual compute (v2, v1) took: " << duration.count() << " microseconds.\n";

//     std::cout<<l2_norm(sx1)<<std::endl;

    // #include <chrono>
    // #include <iostream>

    // // Assuming interpolation_ is already defined and set up
    // auto start = std::chrono::high_resolution_clock::now();
    // interpolation_->applyRestriction(level1, level2, v2, v1);
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "applyRestriction (v2, v1) took: " << duration.count() << " microseconds.\n";

    // start = std::chrono::high_resolution_clock::now();
    // interpolation_->applyRestriction(level1, level2, sv2, sv1);
    // end = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "applyRestriction (sv2, sv1) took: " << duration.count() << " microseconds.\n";

    // start = std::chrono::high_resolution_clock::now();
    // interpolation_->applyFMGInterpolation(level2, level1, v1, v2);
    // end = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "applyFMGInterpolation (v1, v2) took: " << duration.count() << " microseconds.\n";

    // start = std::chrono::high_resolution_clock::now();
    // interpolation_->applyFMGInterpolation(level2, level1, sv1, sv2);
    // end = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "applyFMGInterpolation (sv1, sv2) took: " << duration.count() << " microseconds.\n";


    
}