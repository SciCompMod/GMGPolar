#include "../../../include/Smoother/SmootherTake/smootherTake.h"

#define NODE_BUILD_SMOOTHER_TAKE(i_r, i_theta, grid, DirBC_Interior, inner_boundary_circle_matrix,                     \
                                 circle_tridiagonal_solver, radial_tridiagonal_solver)                                 \
    do {                                                                                                               \
        assert(i_r >= 0 && i_r < grid.nr());                                                                           \
        assert(i_theta >= 0 && i_theta < grid.ntheta());                                                               \
                                                                                                                       \
        const int numberSmootherCircles = grid.numberSmootherCircles();                                                \
        const int lengthSmootherRadial  = grid.lengthSmootherRadial();                                                 \
                                                                                                                       \
        assert(numberSmootherCircles >= 2);                                                                            \
        assert(lengthSmootherRadial >= 3);                                                                             \
                                                                                                                       \
        int row, column;                                                                                               \
        double value;                                                                                                  \
        /* ------------------------------------------ */                                                               \
        /* Node in the interior of the Circle Section */                                                               \
        /* ------------------------------------------ */                                                               \
        if (i_r > 0 && i_r < numberSmootherCircles) {                                                                  \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int left   = grid.index(i_r - 1, i_theta);                                                           \
            const int bottom = grid.index(i_r, i_theta_M1);                                                            \
            const int center = grid.index(i_r, i_theta);                                                               \
            const int top    = grid.index(i_r, i_theta_P1);                                                            \
            const int right  = grid.index(i_r + 1, i_theta);                                                           \
                                                                                                                       \
            auto& matrix           = circle_tridiagonal_solver[i_r];                                                   \
            const int center_index = i_theta;                                                                          \
            const int bottom_index = i_theta_M1;                                                                       \
            const int top_index    = i_theta_P1;                                                                       \
                                                                                                                       \
            /* Center: (Left, Right, Bottom, Top) */                                                                   \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center])                            \
                                                                                                                       \
                     + coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                     coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]));                       \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
                                                                                                                       \
            /* Bottom */                                                                                               \
            row    = center_index;                                                                                     \
            column = bottom_index;                                                                                     \
            value  = -coeff3 * (att[center] + att[bottom]);                                                            \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
                                                                                                                       \
            /* Top */                                                                                                  \
            row    = center_index;                                                                                     \
            column = top_index;                                                                                        \
            value  = -coeff4 * (att[center] + att[top]);                                                               \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Node in the interior of the Radial Section */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r > numberSmootherCircles && i_r < grid.nr() - 2) {                                                 \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int left   = grid.index(i_r - 1, i_theta);                                                           \
            const int bottom = grid.index(i_r, i_theta_M1);                                                            \
            const int center = grid.index(i_r, i_theta);                                                               \
            const int top    = grid.index(i_r, i_theta_P1);                                                            \
            const int right  = grid.index(i_r + 1, i_theta);                                                           \
                                                                                                                       \
            auto& matrix           = radial_tridiagonal_solver[i_theta];                                               \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int left_index   = i_r - numberSmootherCircles - 1;                                                  \
            const int right_index  = i_r - numberSmootherCircles + 1;                                                  \
                                                                                                                       \
            /* Center: (Left, Right, Bottom, Top) */                                                                   \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center])                            \
                                                                                                                       \
                     + coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                     coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]));                       \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
                                                                                                                       \
            /* Left */                                                                                                 \
            row    = center_index;                                                                                     \
            column = left_index;                                                                                       \
            value  = -coeff1 * (arr[center] + arr[left]);                                                              \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
                                                                                                                       \
            /* Right */                                                                                                \
            row    = center_index;                                                                                     \
            column = right_index;                                                                                      \
            value  = -coeff2 * (arr[center] + arr[right]);                                                             \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Circle Section: Node in the inner boundary */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r == 0) {                                                                                           \
            /* ------------------------------------------------ */                                                     \
            /* Case 1: Dirichlet boundary on the inner boundary */                                                     \
            /* ------------------------------------------------ */                                                     \
            if (DirBC_Interior) {                                                                                      \
                auto& matrix           = inner_boundary_circle_matrix;                                                 \
                const int center_index = i_theta;                                                                      \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                const Stencil& CenterStencil = getStencil(i_r);                                                        \
                const int center_nz_index    = getCircleAscIndex(i_r, i_theta);                                        \
                const int nz_index           = center_nz_index + CenterStencil[StencilType::Center];                   \
                matrix.row_index(nz_index)   = center_index + 1;                                                       \
                matrix.col_index(nz_index)   = center_index + 1;                                                       \
                matrix.value(nz_index)       = 1.0;                                                                    \
            }                                                                                                          \
            else {                                                                                                     \
                /* ------------------------------------------------------------- */                                    \
                /* Case 2: Across origin discretization on the interior boundary */                                    \
                /* ------------------------------------------------------------- */                                    \
                /* h1 gets replaced with 2 * R0. */                                                                    \
                /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)). */                          \
                /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */    \
                const double h1 = 2.0 * grid.radius(0);                                                                \
                const double h2 = grid.radialSpacing(i_r);                                                             \
                const double k1 = grid.angularSpacing(i_theta - 1);                                                    \
                const double k2 = grid.angularSpacing(i_theta);                                                        \
                                                                                                                       \
                const double coeff1 = 0.5 * (k1 + k2) / h1;                                                            \
                const double coeff2 = 0.5 * (k1 + k2) / h2;                                                            \
                const double coeff3 = 0.5 * (h1 + h2) / k1;                                                            \
                const double coeff4 = 0.5 * (h1 + h2) / k2;                                                            \
                                                                                                                       \
                const int i_theta_M1           = grid.wrapThetaIndex(i_theta - 1);                                     \
                const int i_theta_P1           = grid.wrapThetaIndex(i_theta + 1);                                     \
                const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);                     \
                                                                                                                       \
                const int left   = grid.index(i_r, i_theta_AcrossOrigin);                                              \
                const int bottom = grid.index(i_r, i_theta_M1);                                                        \
                const int center = grid.index(i_r, i_theta);                                                           \
                const int top    = grid.index(i_r, i_theta_P1);                                                        \
                const int right  = grid.index(i_r + 1, i_theta);                                                       \
                                                                                                                       \
                auto& matrix = inner_boundary_circle_matrix;                                                           \
                                                                                                                       \
                const int center_index = i_theta;                                                                      \
                const int left_index   = i_theta_AcrossOrigin;                                                         \
                const int bottom_index = i_theta_M1;                                                                   \
                const int top_index    = i_theta_P1;                                                                   \
                                                                                                                       \
                const int center_nz_index = getCircleAscIndex(i_r, i_theta);                                           \
                                                                                                                       \
                int nz_index;                                                                                          \
                /* Fill matrix row of (i,j) */                                                                         \
                const Stencil& CenterStencil = getStencil(i_r);                                                        \
                                                                                                                       \
                const double center_value =                                                                            \
                    (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center])                             \
                                                                                                                       \
                     + coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                      \
                     coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]));                        \
                const double left_value   = -coeff1 * (arr[center] + arr[left]);                                       \
                const double bottom_value = -coeff3 * (att[center] + att[bottom]);                                     \
                const double top_value    = -coeff4 * (att[center] + att[top]);                                        \
                                                                                                                       \
                nz_index                   = center_nz_index + CenterStencil[StencilType::Center];                     \
                matrix.row_index(nz_index) = center_index + 1;                                                         \
                matrix.col_index(nz_index) = center_index + 1;                                                         \
                matrix.value(nz_index)     = center_value;                                                             \
                                                                                                                       \
                nz_index                   = center_nz_index + CenterStencil[StencilType::Left];                       \
                matrix.row_index(nz_index) = center_index + 1;                                                         \
                matrix.col_index(nz_index) = left_index + 1;                                                           \
                matrix.value(nz_index)     = left_value;                                                               \
                                                                                                                       \
                nz_index                   = center_nz_index + CenterStencil[StencilType::Bottom];                     \
                matrix.row_index(nz_index) = center_index + 1;                                                         \
                matrix.col_index(nz_index) = bottom_index + 1;                                                         \
                matrix.value(nz_index)     = bottom_value;                                                             \
                                                                                                                       \
                nz_index                   = center_nz_index + CenterStencil[StencilType::Top];                        \
                matrix.row_index(nz_index) = center_index + 1;                                                         \
                matrix.col_index(nz_index) = top_index + 1;                                                            \
                matrix.value(nz_index)     = top_value;                                                                \
            }                                                                                                          \
        }                                                                                                              \
        /* --------------------------------------------- */                                                            \
        /* Radial Section: Node next to circular section */                                                            \
        /* --------------------------------------------- */                                                            \
        else if (i_r == numberSmootherCircles) {                                                                       \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int left   = grid.index(i_r - 1, i_theta);                                                           \
            const int bottom = grid.index(i_r, i_theta_M1);                                                            \
            const int center = grid.index(i_r, i_theta);                                                               \
            const int top    = grid.index(i_r, i_theta_P1);                                                            \
            const int right  = grid.index(i_r + 1, i_theta);                                                           \
                                                                                                                       \
            auto& matrix           = radial_tridiagonal_solver[i_theta];                                               \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int right_index  = i_r - numberSmootherCircles + 1;                                                  \
                                                                                                                       \
            /* Center: (Left, Right, Bottom, Top) */                                                                   \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center])                            \
                                                                                                                       \
                     + coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                     coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]));                       \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
                                                                                                                       \
            /* Right */                                                                                                \
            row    = center_index;                                                                                     \
            column = right_index;                                                                                      \
            value  = -coeff2 * (arr[center] + arr[right]);                                                             \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
        }                                                                                                              \
        /* ------------------------------------------- */                                                              \
        /* Radial Section: Node next to outer boundary */                                                              \
        /* ------------------------------------------- */                                                              \
        else if (i_r == grid.nr() - 2) {                                                                               \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int left   = grid.index(i_r - 1, i_theta);                                                           \
            const int bottom = grid.index(i_r, i_theta_M1);                                                            \
            const int center = grid.index(i_r, i_theta);                                                               \
            const int top    = grid.index(i_r, i_theta_P1);                                                            \
            const int right  = grid.index(i_r + 1, i_theta);                                                           \
                                                                                                                       \
            auto& matrix           = radial_tridiagonal_solver[i_theta];                                               \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int left_index   = i_r - numberSmootherCircles - 1;                                                  \
                                                                                                                       \
            /* Center: (Left, Right, Bottom, Top) */                                                                   \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center])                            \
                                                                                                                       \
                     + coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                     coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]));                       \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
                                                                                                                       \
            /* Left */                                                                                                 \
            row    = center_index;                                                                                     \
            column = left_index;                                                                                       \
            value  = -coeff1 * (arr[center] + arr[left]);                                                              \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Radial Section: Node on the outer boundary */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r == grid.nr() - 1) {                                                                               \
            auto& matrix           = radial_tridiagonal_solver[i_theta];                                               \
            const int center_index = i_r - numberSmootherCircles;                                                      \
                                                                                                                       \
            /* Fill matrix row of (i,j) */                                                                             \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = 1.0;                                                                                              \
            if (row == column)                                                                                         \
                matrix.main_diagonal(row) = value;                                                                     \
            else if (row == column - 1)                                                                                \
                matrix.sub_diagonal(row) = value;                                                                      \
            else if (row == 0 && column == matrix.columns() - 1)                                                       \
                matrix.cyclic_corner_element() = value;                                                                \
        }                                                                                                              \
    } while (0)

void SmootherTake::buildAscCircleSection(const int i_r)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr        = level_cache_.arr();
    const auto& att        = level_cache_.att();
    const auto& art        = level_cache_.art();
    const auto& detDF      = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        // Build Asc at the current node
        NODE_BUILD_SMOOTHER_TAKE(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                                 circle_tridiagonal_solver_, radial_tridiagonal_solver_);
    }
}

void SmootherTake::buildAscRadialSection(const int i_theta)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr        = level_cache_.arr();
    const auto& att        = level_cache_.att();
    const auto& art        = level_cache_.art();
    const auto& detDF      = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        // Build Asc at the current node
        NODE_BUILD_SMOOTHER_TAKE(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                                 circle_tridiagonal_solver_, radial_tridiagonal_solver_);
    }
}

void SmootherTake::buildAscMatrices()
{
    omp_set_num_threads(num_omp_threads_);

    /* -------------------------------------- */
    /* Part 1: Allocate Asc Smoother matrices */
    /* -------------------------------------- */

    const int number_smoother_circles = grid_.numberSmootherCircles();
    const int length_smoother_radial  = grid_.lengthSmootherRadial();

    const int num_circle_nodes = grid_.ntheta();
    circle_tridiagonal_solver_.resize(number_smoother_circles);

    const int num_radial_nodes = length_smoother_radial;
    radial_tridiagonal_solver_.resize(grid_.ntheta());

// Remark: circle_tridiagonal_solver_[0] is unitialized.
// Please use inner_boundary_circle_matrix_ instead!
#pragma omp parallel if (grid_.numberOfNodes() > 10'000)
    {
// ---------------- //
// Circular Section //
#pragma omp for nowait
        for (int circle_Asc_index = 0; circle_Asc_index < number_smoother_circles; circle_Asc_index++) {

            /* Inner boundary circle */
            if (circle_Asc_index == 0) {
                // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
                const int nnz                 = getNonZeroCountCircleAsc(circle_Asc_index);
                inner_boundary_circle_matrix_ = SparseMatrix<double>(num_circle_nodes, num_circle_nodes, nnz);
                inner_boundary_circle_matrix_.is_symmetric(false);
            }

            /* Interior Circle Section */
            else {
                auto& solverMatrix = circle_tridiagonal_solver_[circle_Asc_index];
                solverMatrix       = SymmetricTridiagonalSolver<double>(num_circle_nodes);
                solverMatrix.is_cyclic(true);
            }
        }

// -------------- //
// Radial Section //
#pragma omp for nowait
        for (int radial_Asc_index = 0; radial_Asc_index < grid_.ntheta(); radial_Asc_index++) {
            auto& solverMatrix = radial_tridiagonal_solver_[radial_Asc_index];
            solverMatrix       = SymmetricTridiagonalSolver<double>(num_radial_nodes);
            solverMatrix.is_cyclic(false);
        }
    }

    /* ---------------------------------- */
    /* Part 2: Fill Asc Smoother matrices */
    /* ---------------------------------- */

#pragma omp parallel
    {
#pragma omp for nowait
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            buildAscCircleSection(i_r);
        }

#pragma omp for nowait
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            buildAscRadialSection(i_theta);
        }
    }

    /* ------------------------------------------------------------------ */
    /* Part 3: Convert inner_boundary_circle_matrix to a symmetric matrix */
    /* ------------------------------------------------------------------ */

    SparseMatrix<double> full_matrix = std::move(inner_boundary_circle_matrix_);

    const int nnz           = full_matrix.non_zero_size();
    const int numRows       = full_matrix.rows();
    const int numColumns    = full_matrix.columns();
    const int symmetric_nnz = nnz - (nnz - numRows) / 2;

    inner_boundary_circle_matrix_ = SparseMatrix<double>(numRows, numColumns, symmetric_nnz);
    inner_boundary_circle_matrix_.is_symmetric(true);

    int current_nz = 0;
    for (int nz_index = 0; nz_index < full_matrix.non_zero_size(); nz_index++) {
        int current_row = full_matrix.row_index(nz_index);
        int current_col = full_matrix.col_index(nz_index);
        if (current_row <= current_col) {
            inner_boundary_circle_matrix_.row_index(current_nz) = current_row;
            inner_boundary_circle_matrix_.col_index(current_nz) = current_col;
            inner_boundary_circle_matrix_.value(current_nz)     = std::move(full_matrix.value(nz_index));
            current_nz++;
        }
    }
}