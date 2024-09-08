#include "../../include/Residual/residual.h"

namespace {
    void arr_att_art(const DomainGeometry& domain_geometry,
        const double& r, const double& theta, const double& sin_theta, const double& cos_theta, const double& coeff_alpha, 
        double& arr, double& att, double& art, double& detDF)
    {
        const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta);
        const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta);
        const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta);
        const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta);
        detDF = Jrr * Jtt - Jrt * Jtr;
        arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF);
        att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF);
        art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);
    }
}


void Residual::computeResidualTake0(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const{
    assert(x.size() == grid_.numberOfNodes());
    assert(result.size() == grid_.numberOfNodes());

    result = rhs;

    double scaleAx = -1.0;

    #pragma omp parallel for
    for(int index = 0; index < grid_.numberOfNodes(); index ++){
        double arr, arr_left, arr_right, arr_bottom, arr_top;
        double att, att_left, att_right, att_bottom, att_top;
        double art, art_left, art_right, art_bottom, art_top;
        double detDF, detDF_left, detDF_right, detDF_bottom, detDF_top;

        double coeff_alpha, coeff_beta;
        double sin_theta, cos_theta;

        MultiIndex node = grid_.multiIndex(index);
        Point coords = grid_.polarCoordinates(node);

        coeff_alpha = coeff_alpha_cache_[node[0]];
        coeff_beta = coeff_beta_cache_[node[0]];

        sin_theta = sin(coords[1]); cos_theta = cos(coords[1]);

        arr_att_art(domain_geometry_, coords[0], coords[1], sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF);

        std::array<std::pair<double,double>, space_dimension> neighbor_distance;
        grid_.adjacentNeighborDistances(node, neighbor_distance);

        std::array<std::pair<int,int>, space_dimension> neighbors;
        grid_.adjacentNeighborsOf(node, neighbors);

        if(neighbors[0].second == -1 || (node[0] == 0 && DirBC_Interior_)){
            // Dirichlet Boundary
            result[index] += scaleAx * x[index];
        } else{
            /* Gather arr, art, att values from adjacent neighbors */
            if(neighbors[0].first != -1){
                MultiIndex left_node = grid_.multiIndex(neighbors[0].first);
                Point left_coords = grid_.polarCoordinates(left_node); 
                double coeff_alpha_left = coeff_alpha_cache_[left_node[0]];
                sin_theta = sin(left_coords[1]); cos_theta = cos(left_coords[1]);
                arr_att_art(domain_geometry_, left_coords[0], left_coords[1], sin_theta, cos_theta, coeff_alpha_left, arr_left, att_left, art_left, detDF_left);     
            }else{
                MultiIndex across_origin_node(0, (node[1] + grid_.ntheta() / 2) % grid_.ntheta());
                Point across_origin_coords = grid_.polarCoordinates(across_origin_node);
                double coeff_alpha_left = coeff_alpha_cache_[across_origin_node[0]];
                sin_theta = sin(across_origin_coords[1]); cos_theta = cos(across_origin_coords[1]);
                arr_att_art(domain_geometry_, across_origin_coords[0], across_origin_coords[1], sin_theta, cos_theta, coeff_alpha_left, arr_left, att_left, art_left, detDF_left);  
            }

            // Right
            if(neighbors[0].second != -1){
                MultiIndex right_node = grid_.multiIndex(neighbors[0].second);
                Point right_coords = grid_.polarCoordinates(right_node);
                double coeff_alpha_right = coeff_alpha_cache_[right_node[0]];
                sin_theta = sin(right_coords[1]); cos_theta = cos(right_coords[1]);
                arr_att_art(domain_geometry_, right_coords[0], right_coords[1], sin_theta, cos_theta, coeff_alpha_right, arr_right, att_right, art_right, detDF_right);     
            }
            // Bottom
            if(neighbors[1].first != -1){
                MultiIndex bottom_node = grid_.multiIndex(neighbors[1].first);
                Point bottom_coords = grid_.polarCoordinates(bottom_node);
                double coeff_alpha_bottom = coeff_alpha_cache_[bottom_node[0]];
                sin_theta = sin(bottom_coords[1]); cos_theta = cos(bottom_coords[1]);
                arr_att_art(domain_geometry_, bottom_coords[0], bottom_coords[1], sin_theta, cos_theta, coeff_alpha_bottom, arr_bottom, att_bottom, art_bottom, detDF_bottom);   
            }
            // Top
            if(neighbors[1].second != -1){
                MultiIndex top_node = grid_.multiIndex(neighbors[1].second);
                Point top_coords = grid_.polarCoordinates(top_node);
                double coeff_alpha_top = coeff_alpha_cache_[top_node[0]];
                sin_theta = sin(top_coords[1]); cos_theta = cos(top_coords[1]);
                arr_att_art(domain_geometry_, top_coords[0], top_coords[1], sin_theta, cos_theta, coeff_alpha_top, arr_top, att_top, art_top, detDF_top);   
            }

            double h1 = neighbor_distance[0].first;
            if(node[0] == 0 && !DirBC_Interior_) h1 = 2.0 * grid_.radius(0);
            double h2 = neighbor_distance[0].second;
            double k1 = neighbor_distance[1].first;
            double k2 = neighbor_distance[1].second;

            // beta_{s,t} / f_{s,t} //
            double value = (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF) / 4.0  * x[index];

            // ------ //
            // Center //
            // ------ //
            // Left
            if(neighbors[0].first != -1){
                value += 0.5 * (k1 + k2) / h1 * (arr + arr_left) * x[index];
            }else{
                MultiIndex across_origin_node(0, (node[1] + grid_.ntheta() / 2) % grid_.ntheta());
                value += 0.5 * (k1 + k2) / (2 * grid_.radius(0)) * (arr + arr_left) * x[index];
            }
            // Right
            if(neighbors[0].second != -1){
                value += 0.5 * (k1 + k2) / h2 * (arr + arr_right) * x[index];
            }
            // Bottom
            if(neighbors[1].first != -1){
                value += 0.5 * (h1 + h2) / k1 * (att + att_bottom) * x[index];
            }
            // Top
            if(neighbors[1].second != -1){
                value += 0.5 * (h1 + h2) / k2 * (att + att_top) * x[index];
            }

            // ---- //
            // Left //
            if(neighbors[0].first != -1){
                value += -0.5 * (k1 + k2) / h1 * (arr + arr_left) * x[neighbors[0].first];
            }else{
                MultiIndex across_origin_node(0, (node[1] + grid_.ntheta() / 2) % grid_.ntheta());
                value += -0.5 * (k1 + k2) / (2 * grid_.radius(0)) * (arr + arr_left) * x[grid_.index(across_origin_node)];
            }
            // ----- //
            // Right //
            if(neighbors[0].second != -1){
                value += -0.5 * (k1 + k2) / h2 * (arr + arr_right) * x[neighbors[0].second];
            }
            // ------ //
            // Bottom //
            if(neighbors[1].first != -1){
                value += -0.5 * (h1 + h2) / k1 * (att + att_bottom) * x[neighbors[1].first];
            }
            // --- //
            // Top //
            if(neighbors[1].second != -1){
                value += -0.5 * (h1 + h2) / k2 * (att + att_top) * x[neighbors[1].second];
            }

            grid_.diagonalNeighborsOf(node, neighbors);

            // Bottom Left
            if(neighbors[0].first != -1){
                value += -0.25 * (art_bottom + art_left) * x[neighbors[0].first];
            }
            // Bottom Right
            if(neighbors[0].second != -1){
                value += 0.25 * (art_bottom + art_right) * x[neighbors[0].second];
            }
            // Top Left
            if(neighbors[1].first != -1){
                value += 0.25 * (art_top + art_left) * x[neighbors[1].first];
            }
            // Top Right
            if(neighbors[1].second != -1){
                value += -0.25 * (art_top + art_right) * x[neighbors[1].second];
            }

            result[index] += scaleAx * value;
        }
    }
}