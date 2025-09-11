#include "../../../include/Smoother/SmootherGive/smootherGive.h"

/* ------------------------------------ */
/* Parallelization Version 2: Task Loop */
/* ------------------------------------ */

void SmootherGive::smoothingTaskLoop(Vector<double> x, const Vector<double> rhs, Vector<double> temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

    if (num_omp_threads_ == 1) {
        smoothingSequential(x, rhs, temp);
    }
    else {
        temp = rhs;

        /* Multi-threaded execution */
        const int num_circle_tasks = grid_.numberSmootherCircles();
        const int num_radial_tasks = grid_.ntheta();

        /* Idea: Scedule first_black_circle_smoother as fast as possible to start the radial section early. */
        int first_circle_black_Asc0, first_circle_black_Asc1, first_circle_black_Asc2;
        int first_black_circle_smoother;

        int circle_black_Asc0, circle_black_Asc1, circle_black_Asc2;
        int black_circle_smoother;
        int circle_white_Asc0, circle_white_Asc1, circle_white_Asc2;
        int white_circle_smoother;

        int radial_black_Asc0, radial_black_Asc1, radial_black_Asc2;
        int black_radial_smoother;
        int radial_white_Asc0, radial_white_Asc1, radial_white_Asc2;
        int white_radial_smoother;

#pragma omp parallel num_threads(num_omp_threads_)
        {
            Vector<double> circle_solver_storage_1("circle_solver_storage_1", grid_.ntheta());
            Vector<double> circle_solver_storage_2("circle_solver_storage_2", grid_.ntheta());
            Vector<double> radial_solver_storage("radial_solver_storage", grid_.lengthSmootherRadial());

#pragma omp single
            {
                /* ---------------------------------- */
                /* ------ CIRCLE BLACK SECTION ------ */
                /* ---------------------------------- */

                /* --------------- */
                /* First few lines */

#pragma omp task depend(out : first_circle_black_Asc0)
                {
#pragma omp taskloop
                    for (int circle_task = 0; circle_task < std::min(6, num_circle_tasks); circle_task += 2) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : first_circle_black_Asc0) depend(out : first_circle_black_Asc1)
                {
#pragma omp taskloop
                    for (int circle_task = -1; circle_task < std::min(7, num_circle_tasks); circle_task += 4) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : first_circle_black_Asc1) depend(out : first_circle_black_Asc2)
                {
#pragma omp taskloop
                    for (int circle_task = 1; circle_task < std::min(5, num_circle_tasks); circle_task += 4) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : first_circle_black_Asc2) depend(out : first_black_circle_smoother)
                {
#pragma omp taskloop
                    for (int circle_task = 0; circle_task < std::min(2, num_circle_tasks); circle_task += 2) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
                    }
                }

                /* -------------- */
                /* Leftover lines */

#pragma omp task depend(out : circle_black_Asc0)
                {
#pragma omp taskloop
                    for (int circle_task = 6; circle_task < num_circle_tasks; circle_task += 2) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
                    }
                }

// We don't need depend(in: first_circle_black_Asc0) since there is enough separation.
#pragma omp task depend(in : circle_black_Asc0) depend(out : circle_black_Asc1)
                {
#pragma omp taskloop
                    for (int circle_task = 7; circle_task < num_circle_tasks; circle_task += 4) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : circle_black_Asc1, first_circle_black_Asc1) depend(out : circle_black_Asc2)
                {
#pragma omp taskloop
                    for (int circle_task = 5; circle_task < num_circle_tasks; circle_task += 4) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : circle_black_Asc2, first_circle_black_Asc2) depend(out : black_circle_smoother)
                {
#pragma omp taskloop
                    for (int circle_task = 2; circle_task < num_circle_tasks; circle_task += 2) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
                    }
                }

                /* ---------------------------------- */
                /* ------ CIRCLE WHITE SECTION ------ */
                /* ---------------------------------- */

#pragma omp task depend(in : black_circle_smoother, first_black_circle_smoother) depend(out : circle_white_Asc0)
                {
#pragma omp taskloop
                    for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 2) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : circle_white_Asc0) depend(out : circle_white_Asc1)
                {
#pragma omp taskloop
                    for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 4) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : circle_white_Asc1) depend(out : circle_white_Asc2)
                {
#pragma omp taskloop
                    for (int circle_task = 2; circle_task < num_circle_tasks; circle_task += 4) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : circle_white_Asc2)
                {
#pragma omp taskloop
                    for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 2) {
                        int i_r = num_circle_tasks - circle_task - 1;
                        solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
                    }
                }

                /* ---------------------------------- */
                /* ------ RADIAL BLACK SECTION ------ */
                /* ---------------------------------- */

#pragma omp task depend(in : first_black_circle_smoother) depend(out : radial_black_Asc0)
                {
#pragma omp taskloop
                    for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 2) {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : radial_black_Asc0) depend(out : radial_black_Asc1)
                {
#pragma omp taskloop
                    for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 4) {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : radial_black_Asc1) depend(out : radial_black_Asc2)
                {
#pragma omp taskloop
                    for (int radial_task = 3; radial_task < num_radial_tasks; radial_task += 4) {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : radial_black_Asc2) depend(out : black_radial_smoother)
                {
#pragma omp taskloop
                    for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 2) {
                        int i_theta = radial_task;
                        solveRadialSection(i_theta, x, temp, radial_solver_storage);
                    }
                }

                /* ---------------------------------- */
                /* ------ RADIAL White SECTION ------ */
                /* ---------------------------------- */

#pragma omp task depend(in : black_radial_smoother) depend(out : radial_white_Asc0)
                {
#pragma omp taskloop
                    for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 2) {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : radial_white_Asc0) depend(out : radial_white_Asc1)
                {
#pragma omp taskloop
                    for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 4) {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : radial_white_Asc1) depend(out : radial_white_Asc2)
                {
#pragma omp taskloop
                    for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 4) {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
                    }
                }

#pragma omp task depend(in : radial_white_Asc2)
                {
#pragma omp taskloop
                    for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 2) {
                        int i_theta = radial_task;
                        solveRadialSection(i_theta, x, temp, radial_solver_storage);
                    }
                }
            }
        }
    }
}

/* -------------------------------------------- */
/* Parallelization Version 3: Task Dependencies */
/* -------------------------------------------- */

void SmootherGive::smoothingTaskDependencies(Vector<double> x, const Vector<double> rhs, Vector<double> temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

    if (num_omp_threads_ == 1) {
        smoothingSequential(x, rhs, temp);
    }
    else {
        temp = rhs;

        /* Multi-threaded execution */
        const int num_circle_tasks = grid_.numberSmootherCircles();
        const int num_radial_tasks = grid_.ntheta();

        const int additional_radial_tasks = num_radial_tasks % 3;

        assert(num_circle_tasks >= 2);
        assert(num_radial_tasks >= 3 && num_radial_tasks % 2 == 0);

        /* Make sure to deallocate at the end */
        const int shift           = 3; // Additional space to ensure safe access
        int* asc_ortho_circle_dep = new int[num_circle_tasks + 2 * shift]; // -1 is an additional asc_artho task!
        int* smoother_circle_dep  = new int[num_circle_tasks + 2 * shift];

        int* asc_ortho_radial_dep = new int[num_radial_tasks];
        int* smoother_radial_dep  = new int[num_radial_tasks];

#pragma omp parallel num_threads(num_omp_threads_)
        {
            Vector<double> circle_solver_storage_1("circle_solver_storage_1", grid_.ntheta());
            Vector<double> circle_solver_storage_2("circle_solver_storage_2", grid_.ntheta());
            Vector<double> radial_solver_storage("radial_solver_storage", grid_.lengthSmootherRadial());

#pragma omp single
            {
                /* ---------------------------- */
                /* ------ CIRCLE SECTION ------ */

                /* ---------------------------- */
                /* Asc ortho Black Circle Tasks */
                /* ---------------------------- */

                /* Inside Black Section */
                for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 2) {
#pragma omp task depend(out : asc_ortho_circle_dep[circle_task + shift])
                    {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
                    }
                }
                /* Outside Black Section (Part 1)*/
                for (int circle_task = -1; circle_task < num_circle_tasks; circle_task += 4) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_circle_dep[circle_task + shift])                                                   \
    depend(in                                                                                                          \
           : asc_ortho_circle_dep[circle_task - 1 + shift], asc_ortho_circle_dep[circle_task + 1 + shift])
                    {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
                    }
                }
                /* Outside Black Section (Part 2)*/
                for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 4) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_circle_dep[circle_task + shift])                                                   \
    depend(in                                                                                                          \
           : asc_ortho_circle_dep[circle_task - 2 + shift], asc_ortho_circle_dep[circle_task + 2 + shift])
                    {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
                    }
                }

                /* Black Circle Smoother */
                for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 2) {
#pragma omp task depend(out                                                                                            \
                        : smoother_circle_dep[circle_task + shift])                                                    \
    depend(in                                                                                                          \
           : asc_ortho_circle_dep[circle_task - 1 + shift], asc_ortho_circle_dep[circle_task + 1 + shift])
                    {
                        int i_r = num_circle_tasks - circle_task - 1;
                        solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
                    }
                }

                /* ---------------------------- */
                /* Asc ortho White Circle Tasks */
                /* ---------------------------- */

                /* Inside White Section */
                for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 2) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_circle_dep[circle_task + shift])                                                   \
    depend(in                                                                                                          \
           : smoother_circle_dep[circle_task - 1 + shift], smoother_circle_dep[circle_task + 1 + shift])
                    {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
                    }
                }
                /* Outside White Section (Part 1)*/
                for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 4) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_circle_dep[circle_task + shift])                                                   \
    depend(in                                                                                                          \
           : asc_ortho_circle_dep[circle_task - 1 + shift], asc_ortho_circle_dep[circle_task + 1 + shift])
                    {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
                    }
                }
                /* Outside White Section (Part 2)*/
                for (int circle_task = 2; circle_task < num_circle_tasks; circle_task += 4) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_circle_dep[circle_task + shift])                                                   \
    depend(in                                                                                                          \
           : asc_ortho_circle_dep[circle_task - 2 + shift], asc_ortho_circle_dep[circle_task + 2 + shift])
                    {
                        int i_r = num_circle_tasks - circle_task - 1;
                        applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
                    }
                }

                /* White Circle Smoother */
                for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 2) {
#pragma omp task depend(out                                                                                            \
                        : smoother_circle_dep[circle_task + shift])                                                    \
    depend(in                                                                                                          \
           : asc_ortho_circle_dep[circle_task - 1 + shift], asc_ortho_circle_dep[circle_task + 1 + shift])
                    {
                        int i_r = num_circle_tasks - circle_task - 1;
                        solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
                    }
                }

                /* ---------------------------- */
                /* ------ RADIAL SECTION ------ */

                /* ---------------------------- */
                /* Asc ortho Black Radial Tasks */
                /* ---------------------------- */

                /* Inside Black Section */
                for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 2) {
#pragma omp task depend(out : asc_ortho_radial_dep[radial_task]) depend(in : smoother_circle_dep[0 + shift])
                    {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
                    }
                }
                /* Outside Black Section (Part 1) */
                for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 4) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_radial_dep[radial_task])                                                           \
    depend(in                                                                                                          \
           : asc_ortho_radial_dep[(radial_task - 1 + num_radial_tasks) % num_radial_tasks],                            \
             asc_ortho_radial_dep[(radial_task + 1 + num_radial_tasks) % num_radial_tasks])
                    {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
                    }
                }
                /* Outside Black Section (Part 2) */
                for (int radial_task = 3; radial_task < num_radial_tasks; radial_task += 4) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_radial_dep[radial_task])                                                           \
    depend(in                                                                                                          \
           : asc_ortho_radial_dep[(radial_task - 2 + num_radial_tasks) % num_radial_tasks],                            \
             asc_ortho_radial_dep[(radial_task + 2 + num_radial_tasks) % num_radial_tasks])
                    {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
                    }
                }

                /* Black Radial Smoother */
                for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 2) {
#pragma omp task depend(out                                                                                            \
                        : smoother_radial_dep[radial_task])                                                            \
    depend(in                                                                                                          \
           : asc_ortho_radial_dep[(radial_task - 1 + num_radial_tasks) % num_radial_tasks],                            \
             asc_ortho_radial_dep[(radial_task + 0 + num_radial_tasks) % num_radial_tasks],                            \
             asc_ortho_radial_dep[(radial_task + 1 + num_radial_tasks) % num_radial_tasks])
                    {
                        int i_theta = radial_task;
                        solveRadialSection(i_theta, x, temp, radial_solver_storage);
                    }
                }

                /* ---------------------------- */
                /* Asc ortho White Circle Tasks */
                /* ---------------------------- */

                /* Inside White Section */
                for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 2) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_radial_dep[radial_task])                                                           \
    depend(in                                                                                                          \
           : smoother_radial_dep[(radial_task - 1 + num_radial_tasks) % num_radial_tasks],                             \
             smoother_radial_dep[(radial_task + 0 + num_radial_tasks) % num_radial_tasks],                             \
             smoother_radial_dep[(radial_task + 1 + num_radial_tasks) % num_radial_tasks])
                    {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
                    }
                }
                /* Outside White Section (Part 1) */
                for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 4) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_radial_dep[radial_task])                                                           \
    depend(in                                                                                                          \
           : asc_ortho_radial_dep[(radial_task - 1 + num_radial_tasks) % num_radial_tasks],                            \
             asc_ortho_radial_dep[(radial_task + 1 + num_radial_tasks) % num_radial_tasks])
                    {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
                    }
                }
                /* Outside White Section (Part 2) */
                for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 4) {
#pragma omp task depend(out                                                                                            \
                        : asc_ortho_radial_dep[radial_task])                                                           \
    depend(in                                                                                                          \
           : asc_ortho_radial_dep[(radial_task - 2 + num_radial_tasks) % num_radial_tasks],                            \
             asc_ortho_radial_dep[(radial_task + 2 + num_radial_tasks) % num_radial_tasks])
                    {
                        int i_theta = radial_task;
                        applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
                    }
                }

                /* White Radial Smoother */
                for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 2) {
#pragma omp task depend(out                                                                                            \
                        : smoother_radial_dep[radial_task])                                                            \
    depend(in                                                                                                          \
           : asc_ortho_radial_dep[(radial_task - 1 + num_radial_tasks) % num_radial_tasks],                            \
             asc_ortho_radial_dep[(radial_task + 0 + num_radial_tasks) % num_radial_tasks],                            \
             asc_ortho_radial_dep[(radial_task + 1 + num_radial_tasks) % num_radial_tasks])
                    {
                        int i_theta = radial_task;
                        solveRadialSection(i_theta, x, temp, radial_solver_storage);
                    }
                }
            }
        }
        delete[] asc_ortho_circle_dep;
        delete[] asc_ortho_radial_dep;
        delete[] smoother_circle_dep;
        delete[] smoother_radial_dep;
    }
}
