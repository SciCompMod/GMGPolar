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
    const u_D_Functor& u_D,
    const u_D_Interior_Functor& u_D_Interior
) {
    alpha_ = std::make_shared<alpha_Functor>(alpha);
    beta_ = std::make_shared<beta_Functor>(beta);
    rhs_f_ = std::make_shared<rhs_f_Functor>(rhs_f);
    u_D_ = std::make_shared<u_D_Functor>(u_D);
    u_D_Interior_ = std::make_shared<u_D_Interior_Functor>(u_D_Interior);
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
        // auto [matrixA, rhs] = levels_[current_level].build_system();
        auto [matrixA, rhs] = levels_[current_level].build_symmetric_system();
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Build System time: " << duration_ms.count() << " milliseconds" << std::endl;


        // Vector<double> x(matrixA.columns());
        // assign(x, 3.14);
        // x[17] = 900;
        // x[3] = -199;

        // Vector<double> resultBuild(matrixA.rows());
        // assign(resultBuild, 0.0);
        // multiply(resultBuild, matrixA, x);

        // Vector<double> resultApply(matrixA.rows());
        // assign(resultApply, 0.0);
        // levels_[current_level].applyA(resultApply, x, 1.0);

    
        // for (int i = 0; i < resultApply.size(); i++)
        // {
        //     std::cout<<resultApply[i]<<", "<<resultBuild[i]<<std::endl;
        // }

        const int numThreads = omp_get_max_threads();
        omp_set_num_threads(numThreads);

        std::cout<<numThreads<<std::endl;

        DMUMPS_STRUC_C mumps;
        mumps.job = -1;
        mumps.par = 1;
        mumps.sym = (matrixA.is_symmetric() ? 1 : 0);
        mumps.comm_fortran = -987654;
        dmumps_c(&mumps);

        mumps.ICNTL(1) =  0; // Output stream for error messages.
        mumps.ICNTL(2) =  0; // Output stream for diagnostic printing and statistics local to each MPI process.
        mumps.ICNTL(3) =  0; // Output stream for global information, collected on the host
        mumps.ICNTL(4) =  0; // Level of printing for error, warning, and diagnostic messages.
        mumps.ICNTL(5) =  0; // Controls the matrix input format
        mumps.ICNTL(6)  = 7; // Permutes the matrix to a zero-free diagonal and/or scale the matrix 
        mumps.ICNTL(7) =  5; // Computes a symmetric permutation (ordering) to determine the pivot order to be used for the 
        //                      factorization in case of sequential analysis
        mumps.ICNTL(8) = 77; // Describes the scaling strategy
        mumps.ICNTL(9) =  1; // Computes the solution using A or A^T
        mumps.ICNTL(10) = 0; // Applies the iterative refinement to the computed solution
        mumps.ICNTL(11) = 0; // Computes statistics related to an error analysis of the linear system solved
        mumps.ICNTL(12) = 0; // Defines an ordering strategy for symmetric matrices and is used
        mumps.ICNTL(13) = 0; // Controls the parallelism of the root node
        mumps.ICNTL(14) =    // Controls the percentage increase in the estimated working space
            (matrixA.is_symmetric() ? 5 : 20); 
        mumps.ICNTL(15) = 0; // Exploits compression of the input matrix resulting from a block format
        mumps.ICNTL(16) = 0; // Controls the setting of the number of OpenMP threads
        // ICNTL(17) Doesn't exist
        mumps.ICNTL(18) = 0; // Defines the strategy for the distributed input matrix
        mumps.ICNTL(19) = 0; // Computes the Schur complement matrix
        mumps.ICNTL(20) = 0; // Determines the format (dense, sparse, or distributed) of the right-hand sides
        mumps.ICNTL(21) = 0; // Determines the distribution (centralized or distributed) of the solution vectors.
        mumps.ICNTL(22) = 0; // Controls the in-core/out-of-core (OOC) factorization and solve.
        mumps.ICNTL(23) = 0; // Corresponds to the maximum size of the working memory in MegaBytes that MUMPS can
        //                      allocate per working process
        mumps.ICNTL(24) = 0; // Controls the detection of “null pivot rows”.
        mumps.ICNTL(25) = 0; // Allows the computation of a solution of a deficient matrix and also of a null space basis
        mumps.ICNTL(26) = 0; // Drives the solution phase if a Schur complement matrix has been computed
        mumps.ICNTL(27) = -32; // Controls the blocking size for multiple right-hand sides.
        mumps.ICNTL(28) = 0; // Determines whether a sequential or parallel computation of the ordering is performed
        mumps.ICNTL(29) = 0; // Defines the parallel ordering tool (when ICNTL(28)=1) to be used to compute the fill-in reducing permutation.
        mumps.ICNTL(30) = 0; // Computes a user-specified set of entries in the inverse A^−1 of the original matrix
        mumps.ICNTL(31) = 0; // Indicates which factors may be discarded during the factorization.
        mumps.ICNTL(32) = 0; // Performs the forward elimination of the right-hand sides during the factorization
        mumps.ICNTL(33) = 0; // Computes the determinant of the input matrix.
        mumps.ICNTL(34) = 0; // Controls the conservation of the OOC files during JOB= –3
        mumps.ICNTL(35) = 0; // Controls the activation of the BLR feature
        mumps.ICNTL(36) = 0; // Controls the choice of BLR factorization variant
        mumps.ICNTL(37) = 0; // Controls the BLR compression of the contribution blocks
        mumps.ICNTL(38) = 600; // Estimates compression rate of LU factors
        mumps.ICNTL(39) = 500; // Estimates compression rate of contribution blocks
        // ICNTL(40-47) Don't exist
        mumps.ICNTL(48) = 1; // Multithreading with tree parallelism
        mumps.ICNTL(49) = 0; // Compact workarray id%S at the end of factorization phase
        // ICNTL(50-55) Don't exist
        mumps.ICNTL(56) = 0; // Detects pseudo-singularities during factorization and factorizes the root node with a rankrevealing method
        // ICNTL(57) Doesn't exist
        mumps.ICNTL(58) = 2; // Defines options for symbolic factorization
        // ICNTL(59-60) Don't exist

        mumps.CNTL(1) = -1.0 ; // Relative threshold for numerical pivoting
        mumps.CNTL(2) = -1.0 ; // Stopping criterion for iterative refinement
        mumps.CNTL(3) = 0.0 ; // Determine null pivot rows
        mumps.CNTL(4) = -1.0 ; // Determines the threshold for static pivoting
        mumps.CNTL(5) = 0.0 ; // Defines the fixation for null pivots and is effective only when null pivot row detection is active
        // CNTL(6) Doesn't exist 
        mumps.CNTL(7) = 0.0 ; // Defines the precision of the dropping parameter used during BLR compression
        // CNTL(8-15) Don't exist

        start_time = std::chrono::high_resolution_clock::now();
        mumps.job = 4; 
        assert(matrixA.rows() == matrixA.columns());
        mumps.n = matrixA.rows();
        mumps.nz = matrixA.non_zero_size();
        mumps.irn = matrixA.row_indices_data();
        mumps.jcn = matrixA.column_indices_data();
        mumps.a = matrixA.values_data();

        dmumps_c(&mumps);
        
        end_time = std::chrono::high_resolution_clock::now();
        duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Analysis and Factorization time: " << duration_ms.count() << " milliseconds" << std::endl;

        start_time = std::chrono::high_resolution_clock::now();

        mumps.job = 3; 
        mumps.nrhs = 1;
        mumps.nz_rhs = rhs.size();
        mumps.rhs = rhs.begin();
        mumps.lrhs = rhs.size();

        dmumps_c(&mumps);

        end_time = std::chrono::high_resolution_clock::now();
        duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Solve time: " << duration_ms.count() << " milliseconds" << std::endl;

        if (mumps.info[0] != 0) {
            std::cerr << "Error solving the system: " << mumps.info[0] << std::endl;
        }

        // std::cout<<matrixA<<std::endl;

        // std::cout << "Solution:" << std::endl;
        // for (int i = 0; i < rhs.size(); ++i) {
        //     std::cout << mumps.rhs[i] << std::endl;
        // }


        double error = 0.0;

        double maxi = 0.0;

        const auto& grid = levels_[current_level].grid();

        for (int index = 0; index < rhs.size(); index++)
        {
            MultiIndex node = grid.multiindex(index);
            Point coords = grid.polar_coordinates(node);

            double e = (rhs[index] - (*exact_solution_)(coords[0],coords[1],sin(coords[1]),cos(coords[1]) ));

            // std::cout<<e<<std::endl;

            maxi= std::max(maxi,e);

            error += e*e;

        }
        
        std::cout<<"Error: "<<sqrt(error / rhs.size())<<std::endl;
        std::cout<<"Infity Error: "<<maxi<<std::endl;


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