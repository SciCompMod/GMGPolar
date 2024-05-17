#include "../../include/GMGPolar/gmgpolar.h"

GMGPolar::GMGPolar() : parser_() {
    initializeGrid(); initializeGeometry();
    initializeMultigrid(); initializeGeneral();
    setParameters(0, nullptr);
}

void GMGPolar::setParameters(int argc, char* argv[]) {
    if(argc != 0){
        try {
            parser_.parse_check(argc, argv);
        } catch (const cmdline::cmdline_error &parse_error) {
            std::cerr << "Error: " << parse_error.what() << std::endl;
            std::cerr << "Usage: " << parser_.usage() << std::endl;
        }
    }
    parseGrid(); parseGeometry();
    parseMultigrid(); parseGeneral();
}



void GMGPolar::setup() {
    // Create the finest mesh (level 0)
    PolarGrid finest_grid = createFinestGrid();
    std::cout << "System of size (nr x ntheta) = (" << finest_grid.nr() << " x " << finest_grid.ntheta() << ")\n";
    std::cout << "on the coordinates (r x theta): (" << R0 << ", " << Rmax << ") x (" << 0 << ", " << 2 * M_PI << ")\n";

    // Building the grid on all levels
    numberOflevels_ = numberOfLevels(finest_grid);
    initializeLevels(numberOflevels_, finest_grid);
    std::cout<<"Number of Levels used: "<<numberOflevels_<<std::endl;

    interpolation_ = std::make_unique<Interpolation>();
}


void GMGPolar::solve() {


// Vector<scalar_t> result1(getLevel(0).grid().number_of_nodes());
// Vector<scalar_t> result2(getLevel(0).grid().number_of_nodes());
//     {
//         const int current_level = 0;
//         Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

//         for (int i = 0; i < x.size(); i++)
//         {
//             x[i] = 4 * i;
//         }

//         auto start_time = std::chrono::high_resolution_clock::now();
        
//         levels_[current_level].applyATake0(result1, x);

//         auto end_time = std::chrono::high_resolution_clock::now();

//         // Compute duration in milliseconds
//         auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

//         // Output the duration
//         std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

//         // std::cout<<"Input: \n"<<x<<std::endl;
//         // std::cout<<"Result: \n"<<result<<std::endl;
//     }



//     {
//         const int current_level = 0;
//         Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());
//         Vector<scalar_t> result(getLevel(current_level).grid().number_of_nodes());

//         for (int i = 0; i < x.size(); i++)
//         {
//             x[i] = 4 * i;
//         }

//         auto start_time = std::chrono::high_resolution_clock::now();
        
//         levels_[current_level].applyA(result2, x);

//         auto end_time = std::chrono::high_resolution_clock::now();

//         // Compute duration in milliseconds
//         auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

//         // Output the duration
//         std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

//         // std::cout<<"Input: \n"<<x<<std::endl;
//         std::cout<<"Result: \n"<<result<<std::endl;
//     }

//     for (int i = 0; i < result1.size(); i++)
//     {
//         std::cout<<result1[i]<<", "<<result2[i]<<std::endl;
//     }
    

}


const Level& GMGPolar::getLevel(const int index) const {
    assert(index >= 0 && index < numberOflevels_);
    return levels_[index];
}

void GMGPolar::prolongateToNextLevel(const int current_level, Vector<scalar_t> &result, const Vector<scalar_t> &x) const {
    assert(static_cast<size_t>(current_level) < levels_.size() && 1 <= current_level);
    interpolation_->applyProlongation(levels_[current_level], levels_[current_level-1], result, x);
}

void GMGPolar::restrictToLowerLevel(const int current_level, Vector<scalar_t> &result, const Vector<scalar_t> &x) const {
    assert(static_cast<size_t>(current_level) < levels_.size() - 1 && 0 <= current_level);
    interpolation_->applyRestrictionTake(levels_[current_level], levels_[current_level+1], result, x);
}
    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     applyA(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }



    // {
    //     const int current_level = 1;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     applyA(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }


// }








    // int S1_0 = getLevel(0).grid().numberSmootherCircles();
    // std::cout<< getLevel(0).grid().radius(S1_0-1)<< ", "<< getLevel(0).grid().radius(S1_0) <<std::endl;
    // std::cout<<getLevel(0).grid().smoother_splitting_radius()<<std::endl;

    // int S1_1 = getLevel(1).grid().numberSmootherCircles();
    // std::cout<< getLevel(1).grid().radius(S1_1-1)<< ", "<< getLevel(1).grid().radius(S1_1) <<std::endl;
    // std::cout<<getLevel(1).grid().smoother_splitting_radius()<<std::endl;

    // int S1_2 = getLevel(2).grid().numberSmootherCircles();
    // std::cout<< getLevel(2).grid().radius(S1_2-1)<< ", "<< getLevel(2).grid().radius(S1_2) <<std::endl;
    // std::cout<<getLevel(2).grid().smoother_splitting_radius()<<std::endl;


    // op_ = std::make_unique<Operator>();


    // Restriction






    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTake0(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }


    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelGive(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }




    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelGiveTasks(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }





    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTake0(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }




    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTake(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }




    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTakeTasks(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }


    // std::cout<<"AGAIN"<<std::endl;


    //     {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTake0(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }


    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelGive(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }




    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelGiveTasks(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }





    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTake0(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }




    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTake(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }




    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTakeTasks(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }



//    {
//         const int current_level = 0;
//         Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

//         std::random_device rd;
//         std::mt19937 gen(42);
//         std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

//         // for (scalar_t& val : x) {
//         //     val = dis(gen);
//         // }

//         for (int i = 0; i < x.size(); i++)
//         {
//             x[i] = 4 * i;
//         }

//         // Call the function
//         Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
//         auto start_time = std::chrono::high_resolution_clock::now();
//         restrictToLowerLevelTake(current_level, result, x);

//         // End time
//         auto end_time = std::chrono::high_resolution_clock::now();

//         // Compute duration in milliseconds
//         auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

//         // Output the duration
//         std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;


//         //         std::cout<<"Input: \n"<<x<<std::endl;
//         // std::cout<<"Result: \n"<<result<<std::endl;
//     }



//    {
//         const int current_level = 0;
//         Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

//         std::random_device rd;
//         std::mt19937 gen(42);
//         std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

//         // for (scalar_t& val : x) {
//         //     val = dis(gen);
//         // }

//         for (int i = 0; i < x.size(); i++)
//         {
//             x[i] = 4 * i;
//         }

//         // Call the function
//         Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
//         auto start_time = std::chrono::high_resolution_clock::now();
//         restrictToLowerLevelGive0(current_level, result, x);

//         // End time
//         auto end_time = std::chrono::high_resolution_clock::now();

//         // Compute duration in milliseconds
//         auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

//         // Output the duration
//         std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;


//         //         std::cout<<"Input: \n"<<x<<std::endl;
//         // std::cout<<"Result: \n"<<result<<std::endl;
//     }




    
//    {
//         const int current_level = 0;
//         Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

//         std::random_device rd;
//         std::mt19937 gen(42);
//         std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

//         // for (scalar_t& val : x) {
//         //     val = dis(gen);
//         // }

//         for (int i = 0; i < x.size(); i++)
//         {
//             x[i] = 4 * i;
//         }

//         // Call the function
//         Vector<scalar_t> result(getLevel(1).grid().number_of_nodes());
//         auto start_time = std::chrono::high_resolution_clock::now();
//         restrictToLowerLevelGive(current_level, result, x);

//         // End time
//         auto end_time = std::chrono::high_resolution_clock::now();

//         // Compute duration in milliseconds
//         auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

//         // Output the duration
//         std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;


//         // std::cout<<"Input: \n"<<x<<std::endl;
//         // std::cout<<"Result: \n"<<result<<std::endl;
//     }







    // {
    //     const int current_level = 0;
    //     std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     std::vector<scalar_t> result;
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTake0(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // disp(x, "Input");
    //     // disp(result,"Result");
    // }





    // {
    //     const int current_level = 0;
    //     std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     std::vector<scalar_t> result;
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelTake(current_level, result, x);

    //     std::cout<<"Test"<<std::endl;

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // disp(x, "Input");
    //     // disp(result,"Result");
    // }






    // {
    //     const int current_level = 0;
    //     std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     std::vector<scalar_t> result;
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelGive(current_level, result, x);

    //     std::cout<<"Test"<<std::endl;

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // disp(x, "Input");
    //     // disp(result,"Result");
    // }




    // {
    //     const int current_level = 0;
    //     std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     std::vector<scalar_t> result;
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelGiveAtomic(current_level, result, x);

    //     std::cout<<"Test"<<std::endl;

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // disp(x, "Input");
    //     // disp(result,"Result");
    // }



    // {
    //     const int current_level = 0;
    //     std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     std::vector<scalar_t> result;
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     restrictToLowerLevelGiveTasks(current_level, result, x);

    //     std::cout<<"Test"<<std::endl;

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // disp(x, "Input");
    //     // disp(result,"Result");
    // }



   // Prolongation

    // {
    //     const int current_level = 1;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(0).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     prolongateToNextLevel(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // disp(x, "Input");
    //     // disp(result,"Result");

    //     //         std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }



    // {
    //     const int current_level = 1;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(0).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     prolongateToNextLevel(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // disp(x, "Input");
    //     // disp(result,"Result");

    //     //         std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }



    // {
    //     const int current_level = 1;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     // Call the function
    //     Vector<scalar_t> result(getLevel(0).grid().number_of_nodes());
    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     prolongateToNextLevelFor(current_level, result, x);

    //     // End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // disp(x, "Input");
    //     // disp(result,"Result");

    //     //         std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }




 //  Prolongation

    // {
    //     const int current_level = 1;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     std::random_device rd;
    //     std::mt19937 gen(42);
    //     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    //     // for (scalar_t& val : x) {
    //     //     val = dis(gen);
    //     // }

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //    // Call the function
    //     Vector<scalar_t> result(getLevel(0).grid().number_of_nodes());

    //     auto start_time = std::chrono::high_resolution_clock::now();
    //     prolongateToNextLevel0(current_level, result, x);

    //     //End time
    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     //Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     //Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;

    //     // disp(x, "Input");
    //     // disp(result,"Result");
    // }








// {

//     const int current_level = 0;
//     std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

//     std::random_device rd;
//     std::mt19937 gen(42);
//     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

//     // for (scalar_t& val : x) {
//     //     val = dis(gen);
//     // }

//     for (int i = 0; i < x.size(); i++)
//     {
//         x[i] = 4 * i;
//     }

//     // Call the function
//     std::vector<scalar_t> result(getLevel(1).grid().number_of_nodes(), 0);
//     auto start_time = std::chrono::high_resolution_clock::now();
//     restrictToLowerLevel(current_level, result, x);

//     // End time
//     auto end_time = std::chrono::high_resolution_clock::now();

//     // Compute duration in milliseconds
//     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

//     // Output the duration
//     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

//     disp(x, "Input");
//     disp(result,"Result");



// }





//     {

    // const int current_level = 0;
    // std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

    // std::random_device rd;
    // std::mt19937 gen(42);
    // std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    // // for (scalar_t& val : x) {
    // //     val = dis(gen);
    // // }

    // for (int i = 0; i < x.size(); i++)
    // {
    //     x[i] = 7.0;
    // }
    

    // // Call the function
    // std::vector<scalar_t> result(getLevel(1).grid().number_of_nodes(), 0);
    // auto start_time = std::chrono::high_resolution_clock::now();
    // restrictToLowerLevel(current_level, result, x);

//     // End time
//     auto end_time = std::chrono::high_resolution_clock::now();

//     // Compute duration in milliseconds
//     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

//     // Output the duration
//     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

//     disp(x, "Input");
//     disp(result,"Result");
// }


// {

//     const int current_level = 1;
//     std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

//     std::random_device rd;
//     std::mt19937 gen(42);
//     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

//     // for (scalar_t& val : x) {
//     //     val = dis(gen);
//     // }

//     // for (int i = 0; i < x.size(); i++)
//     // {
//     //     x[i] = 4 * i;
//     // }

//     // Call the function
//     std::vector<scalar_t> result(getLevel(0).grid().number_of_nodes(), 0);
//     auto start_time = std::chrono::high_resolution_clock::now();
//     prolongateToNextLevel2(current_level, result, x);

//     // End time
//     auto end_time = std::chrono::high_resolution_clock::now();

//     // Compute duration in milliseconds
//     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

//     // Output the duration
//     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

//     // disp(x, "Input");
//     // disp(result,"Result");
// }





// {

//     const int current_level = 1;
//     std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

//     std::random_device rd;
//     std::mt19937 gen(42);
//     std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

//     // for (scalar_t& val : x) {
//     //     val = dis(gen);
//     // }

//     // for (int i = 0; i < x.size(); i++)
//     // {
//     //     x[i] = 4 * i;
//     // }

//     // Call the function
//     std::vector<scalar_t> result(getLevel(0).grid().number_of_nodes(), 0);
//     auto start_time = std::chrono::high_resolution_clock::now();
//     prolongateToNextLevel(current_level, result, x);

//     // End time
//     auto end_time = std::chrono::high_resolution_clock::now();

//     // Compute duration in milliseconds
//     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

//     // Output the duration
//     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

//     // disp(x, "Input");
//     // disp(result,"Result");
// }








    // const int current_level = 0;
    // std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

    // std::random_device rd;
    // std::mt19937 gen(42);
    // std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    // for (scalar_t& val : x) {
    //     val = dis(gen);
    // }


    // Call the function
    // std::vector<scalar_t> result(getLevel(1).grid().number_of_nodes(), 0);
    // auto start_time = std::chrono::high_resolution_clock::now();
    // restrictToNextLevel(current_level, result, x);

    // End time
    // auto end_time = std::chrono::high_resolution_clock::now();

    // Compute duration in milliseconds
    // auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    // Output the duration
    // std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    // disp(x, "Input");
    // disp(result,"Result");



    // const int current_level = 1;
    // std::vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes(), 0);

    // std::random_device rd;
    // std::mt19937 gen(42);
    // std::uniform_real_distribution<scalar_t> dis(0.0, 1.0); // Range [0.0, 1.0), you can adjust it

    // for (scalar_t& val : x) {
    //     val = dis(gen);
    // }


    // // Call the function
    // std::vector<scalar_t> result(getLevel(0).grid().number_of_nodes(), 0);
    // auto start_time = std::chrono::high_resolution_clock::now();
    // prolongateToNextLevel(current_level, result, x);

    // // End time
    // auto end_time = std::chrono::high_resolution_clock::now();

    // // Compute duration in milliseconds
    // auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    // // Output the duration
    // std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    // // disp(x, "Input");
    // // disp(result,"Result");

// }