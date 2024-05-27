#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::setRadialRefinement(double r_jump){
    r_jump_ = r_jump;
}

void GMGPolar::setGeometry(
    const dFx_dr_Functor& dFx_dr, 
    const dFy_dr_Functor& dFy_dr, 
    const dFx_dt_Functor& dFx_dt, 
    const dFy_dt_Functor& dFy_dt
) {
    dFx_dr_ = std::make_shared<dFx_dr_Functor>(dFx_dr);
    dFy_dr_ = std::make_shared<dFy_dr_Functor>(dFy_dr);
    dFx_dt_ = std::make_shared<dFx_dt_Functor>(dFx_dt);
    dFy_dt_ = std::make_shared<dFy_dt_Functor>(dFy_dt);
}

void GMGPolar::setParameters(
    const alpha_Functor& alpha, 
    const beta_Functor& beta, 
    const rhs_f_Functor& rhs_f, 
    const u_D_Functor& u_D
) {
    alpha_ = std::make_shared<alpha_Functor>(alpha);
    beta_ = std::make_shared<beta_Functor>(beta);
    rhs_f_ = std::make_shared<rhs_f_Functor>(rhs_f);
    u_D_ = std::make_shared<u_D_Functor>(u_D);
}

void GMGPolar::setSystemParameters(const exact_solution_Functor& exact_solution) {
    exact_solution_ = std::make_shared<exact_solution_Functor>(exact_solution);
}

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



    {
        const int current_level = 0;

        auto start_time = std::chrono::high_resolution_clock::now();

        const auto [matrixA, rhs] = levels_[current_level].build_system();

        auto end_time = std::chrono::high_resolution_clock::now();

        // Compute duration in milliseconds
        auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        // Output the duration
        std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

        // std::cout<<"Input: \n"<<x<<std::endl;
        // std::cout<<"Result: \n"<<result<<std::endl;
    }



    // {
    // const int numThreads = omp_get_max_threads();

    // const int NUM_LEVELS = 1;
    // const int NUM_REPEATS = 1;

    // // Store timings for each level
    // std::vector<std::vector<int>> all_timings(NUM_LEVELS);

    // for (int current_level = 0; current_level < NUM_LEVELS; ++current_level) {
    //     std::vector<int> timings;

    //     Vector<scalar_t> result(getLevel(current_level).grid().number_of_nodes());
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     #pragma omp parallel for
    //     for (int i = 0; i < x.size(); i++) { 
    //         x[i] = 4 * i;
    //     }

    //     for (int i = numThreads; i <= numThreads; i += 1) {
    //         omp_set_num_threads(i);

    //         assign(result, 0.0);

    //         std::vector<int> duration_repeats(NUM_REPEATS);

    //         for (int j = 0; j < NUM_REPEATS; ++j) {
                
    //             auto start_time = std::chrono::high_resolution_clock::now();
                
    //             levels_[current_level].applyA(result, x, 1.0);
    //             // levels_[current_level].applyATasks(result, x, 1.0);
    //             // levels_[current_level].applyAMutex(result, x, 1.0);
    //             // levels_[current_level].applyATake0(result, x, 1.0);

    //             auto end_time = std::chrono::high_resolution_clock::now();

    //             int duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    //             duration_repeats[j] = duration;
    //         }

    //         // Compute average duration
    //         int avg_duration = std::accumulate(duration_repeats.begin(), duration_repeats.end(), 0) / NUM_REPEATS;
    //         timings.push_back(avg_duration);

    //         // Print timings to console
    //         std::cout << "Level: " << current_level << " Threads: " << i << " Average Duration: " << avg_duration << " microseconds" << std::endl;
    //     }

    //     all_timings[current_level] = timings;
    // }

    // // Save timings to a file
    // std::ofstream timings_file("GiveMutexTasking2.txt");
    // for (int level = 0; level < NUM_LEVELS; ++level) {
    //     timings_file << "Level " << level << "\n";
    //     for (const auto& timing : all_timings[level]) {
    //         timings_file << timing << "\n";
    //     }
    //     timings_file << "\n";  // Separate levels by a blank line
    // }
    // timings_file.close();

  
    // }


    // Vector<scalar_t> result1(getLevel(0).grid().number_of_nodes());
    // Vector<scalar_t> result2(getLevel(0).grid().number_of_nodes());
    // Vector<scalar_t> result3(getLevel(0).grid().number_of_nodes());
    // Vector<scalar_t> result4(getLevel(0).grid().number_of_nodes());

    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     auto start_time = std::chrono::high_resolution_clock::now();

    //     assign(result1, 0.0);
    //     double factor = 1.0;
        
    //     levels_[current_level].applyATake0(result1, x, factor);

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
    //     // Vector<scalar_t> result(getLevel(current_level).grid().number_of_nodes());

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     auto start_time = std::chrono::high_resolution_clock::now();

    //     assign(result2, 0.0);

    //     double factor = 1.0;
        
    //     levels_[current_level].applyA(result2, x, factor);

    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result2<<std::endl;
    // }


    // {
    //     const int current_level = 0;
    //     Vector<scalar_t> x(getLevel(current_level).grid().number_of_nodes());

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     assign(result3, 0.0);

    //     double factor = 1.0;

    //     auto start_time = std::chrono::high_resolution_clock::now();
        
    //     levels_[current_level].applyATasks(result3, x, factor);

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

    //     for (int i = 0; i < x.size(); i++)
    //     {
    //         x[i] = 4 * i;
    //     }

    //     assign(result4, 0.0);

    //     double factor = 1.0;
        

    //     auto start_time = std::chrono::high_resolution_clock::now();
        
    //     levels_[current_level].applyAMutex(result4, x, factor);

    //     auto end_time = std::chrono::high_resolution_clock::now();

    //     // Compute duration in milliseconds
    //     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    //     // Output the duration
    //     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

    //     // std::cout<<"Input: \n"<<x<<std::endl;
    //     // std::cout<<"Result: \n"<<result<<std::endl;
    // }


    // for (int i = 0; i < result2.size(); i++)
    // {
    //    //  std::cout<<result1[i]<<", "<<result2[i]<<", "<<result4[i]<<std::endl;
    //     std::cout<<result1[i]<<", "<<result2[i]<<", "<<result3[i]<<", "<<result4[i]<<std::endl;
    // }
    

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