#include "mockgrid.h"
void create_grid(gmgpolar& test_p)
{
    std::default_random_engine gen(time(0));
    std::uniform_real_distribution<double> dis(gyro::dcntl[Param::R0], gyro::dcntl[Param::R]);
    std::uniform_real_distribution<double> theta_distribution(0, 2 * PI);
    level* new_level = new level(0);
    new_level->nr    = pow(2, gyro::icntl[Param::nr_exp]);
    new_level->r     = std::vector<double>(new_level->nr + 1);

    std::vector<double> rands(new_level->nr - 1);
    for (int j = 0; j < new_level->nr - 1; ++j) {
        rands[j] = dis(gen);
    }
    std::sort(rands.begin(), rands.end());
    new_level->r[0] = gyro::dcntl[Param::R0];
    for (int i = 1; i < new_level->nr; i++) {
        new_level->r[i] = rands[i - 1]; //random grid on coarser level
    }
    new_level->r[new_level->nr] = gyro::dcntl[Param::R];
    new_level->nr++;
    int ntmp          = pow(2, ceil(log2(new_level->nr)));
    new_level->ntheta = gyro::icntl[Param::periodic] ? ntmp : ntmp + 1;

    new_level->theta = std::vector<double>(new_level->ntheta);

    std::vector<double> randst(ntmp - 1);
    for (int k = 0; k < ntmp - 1; ++k) {
        randst[k] = theta_distribution(gen);
    }
    std::sort(randst.begin(), randst.end());
    new_level->theta[0] = 0;
    for (int i = 1; i < ntmp; i++) {
        new_level->theta[i] = randst[i - 1];
    }
    if (!gyro::icntl[Param::periodic]) {
        new_level->theta[ntmp] = 2 * PI;
    }

    new_level->ntheta_int = gyro::icntl[Param::periodic] ? new_level->ntheta : new_level->ntheta - 1;
    new_level->nr_int     = new_level->nr - 1;

    new_level->thetaplus = std::vector<double>(new_level->ntheta_int);
    for (int k = 0; k < new_level->ntheta_int - 1; k++) {
        new_level->thetaplus[k] = new_level->theta[k + 1] - new_level->theta[k];
    }
    new_level->thetaplus[new_level->ntheta_int - 1] = 2 * PI - new_level->theta[new_level->ntheta_int - 1];

    new_level->hplus = std::vector<double>(new_level->nr_int);
    for (int k = 0; k < new_level->nr_int; k++) {
        new_level->hplus[k] = new_level->r[k + 1] - new_level->r[k];
    }

    level* coarser_level  = new level(1);
    coarser_level->nr     = pow(2, gyro::icntl[Param::nr_exp] - 1) + 1;
    ntmp                  = pow(2, ceil(log2(coarser_level->nr)));
    coarser_level->ntheta = gyro::icntl[Param::periodic] ? ntmp : ntmp + 1;

    coarser_level->ntheta_int = gyro::icntl[Param::periodic] ? coarser_level->ntheta : coarser_level->ntheta - 1;
    coarser_level->nr_int     = coarser_level->nr - 1;
    test_p.v_level.push_back(new_level);
    test_p.v_level.push_back(coarser_level);
    test_p.levels_orig = 2;
}