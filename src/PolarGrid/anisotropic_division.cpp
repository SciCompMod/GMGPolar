#include "../../include/PolarGrid/polargrid.h"

void PolarGrid::RadialAnisotropicDivision(std::vector<double>& r_temp, 
    const double R0, const double R, const int nr_exp, 
    const int anisotropic_factor, const alpha_coeff alpha) const 
{
    // 1) uniform division with nr=2^dummy_lognr - 2^aniso
    // 2) remaining nodes are added by refining the part centered around 2/3 of r
    std::set<double, std::greater<double>>::iterator itr, itr_p1;
    // very ugly anisotropy hack.... dividing recursively smaller and smaller number of cells

    /* uniform division of r in 2^nr_exp - 2^aniso */
    int dummy_lognr  = nr_exp;
    int n_elems_equi = pow(2, dummy_lognr) - pow(2, anisotropic_factor);
    if (anisotropic_factor < 0 || n_elems_equi <= 0) {
        throw std::runtime_error("Please choose anisotropy factor a such that 2^fac_ani < 2^nr_exp.\n");
    }

    if ((anisotropic_factor % 2) == 1) // odd number of elements on an open circular disk is desired because of coarsening
        n_elems_equi++;
    double uniform_distance = (R - R0) / n_elems_equi;
    int nr = n_elems_equi + 1;
    std::vector<double> r_temp2 = std::vector<double>(nr);
    for (int i = 0; i < nr - 1; i++)
        r_temp2[i] = R0 + i * uniform_distance;
    r_temp2[nr - 1] = R;

    /* refine around 2/3 of r */
    int n_elems_refined = pow(2, anisotropic_factor);

    // edge
    int se;

    // allow refining of the grid at r_jump, the center point of the 
    // drop of the diffusion coefficient alpha.
    double r_jump;
    if (alpha == SONNENDRUCKER) {
        // The center of the coefficient jump lies at 0.6888697651782026
        // for backward stability with previous runs and the Matlab code, 
        // we use 0.66 though.
        r_jump = 0.66;
    } else if (alpha == ZONI) {
        r_jump = 0.4837;
    } else if (alpha == ZONI_SHIFTED) {
        // Choose center point of descent.
        // a) - ln(0.5 * (alpha(0) - alpha(Rmax))):
        //    - ln(0.5 * (np.exp(-np.tanh(-14)) - np.exp(-np.tanh(6)))) = 0.16143743821247852
        // b) r_center = Rmax * (np.arctanh(0.16143743821247852) + 14) / 20 = 0.7081431124450334 Rmax
        r_jump = 0.7081;
    } else if (alpha == POISSON) {
        r_jump = 0.5; // There is no jump for Poisson so this is an arbitrary choice
    } else {
        throw std::runtime_error("Unknown alpha coeff");
    }

    // Added by Allan Kuhn to fix a memory error
    if (floor(nr * r_jump) > nr - (n_elems_refined / 2)) {
        int new_aniso   = log2(nr - floor(nr * r_jump)) + 1;
        n_elems_refined = pow(2, new_aniso);
    } 

    se = floor(nr * r_jump) - n_elems_refined / 2;
    int ee = se + n_elems_refined;
    // takeout
    int st = ceil((double)n_elems_refined / 4.0 + 1) - 1;
    int et = floor(3 * ((double)n_elems_refined / 4.0));

    std::set<double> r_set;
    std::set<double> r_set_p1;
    int count = 0;
    for (int i = 0; i < n_elems_refined; i++) {
        r_set_p1.insert(r_temp2[se + i]);
        count++;
    }
    double half = uniform_distance / 2.0;
    for (int k = 0; k < anisotropic_factor; k++) {
        std::set<double> r_set_p1_tmp;
        itr_p1     = r_set_p1.begin();
        int r_size = count;
        count      = 0;
        for (int i = 0; i < r_size - 1; i++) {
            r_set.insert((*itr_p1) + half);
            if (k < anisotropic_factor - 1 && i >= st && i < et) {
                r_set_p1_tmp.insert(*(itr_p1));
                r_set_p1_tmp.insert(*(itr_p1) + half);
                count += 2;
            }
            itr_p1++;
        }
        r_set_p1 = r_set_p1_tmp;
        half *= 0.5;
    }

    // such that the total size is 8*x+1 (or we do not refine)
    nr        = nr + r_set.size();
    int shift = 0;
    shift     = std::min(nr % 8 - 1, (int)r_set.size());
    itr       = r_set.begin();
    std::advance(itr, shift);
    r_set.erase(r_set.begin(), itr);
    for (int i = 0; i < n_elems_refined; i++)
        r_set.insert(r_temp2[se + i]);

    // group all in r_tmp
    nr = n_elems_equi - n_elems_refined + r_set.size() + 1;

    r_temp.resize(nr);

    for (int i = 0; i < se; i++)
        r_temp[i] = r_temp2[i];
    itr = r_set.begin();
    for (int i = 0; i < (int)r_set.size(); i++) {
        r_temp[se + i] = *itr;
        itr++;
    }
    for (int i = 0; i < n_elems_equi - ee + 1; i++)
        r_temp[se + r_set.size() + i] = r_temp2[ee + i];
}