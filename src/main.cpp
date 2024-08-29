#include <iostream>
#include <chrono>
#include <memory>
#include <random>

#include "../include/GMGPolar/gmgpolar.h"

#include <chrono>
#include <thread>

int main(int argc, char* argv[]){
    #ifdef NDEBUG
        std::cout << "Build Type: Release\n"<<std::endl;
    #else
        std::cout << "Build Type: Debug\n"<<std::endl;
    #endif

    GMGPolar solver;

    // Configure solver parameters from command-line arguments
    solver.setParameters(argc, argv);
    
    // Setup and solve the problem
    solver.setup();
    solver.solve();

    // Recover solution
    Vector<double>& solution = solver.solution();
    const PolarGrid& grid = solver.grid();

    // Print Timings
    solver.printTimings();



    // int A_done = 0;
    // int C_done = 0;

    // auto start = std::chrono::high_resolution_clock::now();

    // #pragma omp parallel
    // {
    //     #pragma omp single
    //     {
    //         // Task for loop A
    //         #pragma omp task depend(out: A_done)
    //         {
    //             #pragma omp taskloop
    //             for (int i = 0; i < 3; i++) {
    //                 std::this_thread::sleep_for(std::chrono::seconds(1));
    //                 std::cout << "A: " << i << std::endl;
    //             }
    //         }

    //         // Task for loop B dependent on A_done
    //         #pragma omp task depend(in: A_done)
    //         {
    //             #pragma omp taskloop
    //             for (int i = 0; i < 3; i++) {
    //                 std::this_thread::sleep_for(std::chrono::seconds(2));
    //                 std::cout << "B: " << i << std::endl;
    //             }
    //         }

    //         // Task for loop C
    //         #pragma omp task depend(out: C_done)
    //         {
    //             #pragma omp taskloop
    //             for (int i = 0; i < 3; i++) {
    //                 std::this_thread::sleep_for(std::chrono::seconds(4));
    //                 std::cout << "C: " << i << std::endl;
    //             }
    //         }

    //         // Task for loop D dependent on C_done
    //         #pragma omp task depend(in: C_done)
    //         {
    //             #pragma omp taskloop
    //             for (int i = 0; i < 3; i++) {
    //                 std::this_thread::sleep_for(std::chrono::seconds(6));
    //                 std::cout << "D: " << i << std::endl;
    //             }
    //         }
    //     }
    // }

    

    // auto end = std::chrono::high_resolution_clock::now();


    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


    // std::cout << "Operation took " << duration << " milliseconds." << std::endl;







    return 0;
}
