#ifndef LENNARD_JONES_POTENTIAL_H
#define LENNARD_JONES_POTENTIAL_H

#include <string>
#include <vector>
#include <cmath>
#include <iostream>



class LJPotential
{
public:
    LJPotential();

    // Methods
    void ReadLJParameters(const std::string& filename);
    float CalcLJPotential(
        std::vector<std::vector<float>> pdb2gmx_coord_input,
        std::vector<std::vector<float>>dist_matrix_input,
        std::vector<std::string> pdb2gmx_atomlist_input
     );
    std::vector<std::vector<float>> PairwiseCalcLJPotential(
        std::vector<std::vector<float>> pdb2gmx_coord_input,
        std::vector<std::string> pdb2gmx_atomlist_input);

    float PBCCalcLJPotential(
                        std::vector<std::vector<float>> position_inp,
                        float sigma,
                        float epsilon,
                        bool pbc,
                        float box_length);

    std::vector<float> CalcForceLJPotential(
        std::vector<std::vector<float>> pdb2gmx_coord_input,
        std::vector<std::string> pdb2gmx_atomlist_input,
        bool pbc,
        float box_length_input);

    std::vector<float> PBCForceLJPotential(
                        std::vector<std::vector<float>> position_inp,
                        float sigma,
                        float epsilon,
                        bool pbc,
                        float box_length);

    //Variables

    std::vector<std::string> lj_atompair;
    std::vector<float> lj_sigma;
    std::vector<float> lj_epsilon;
    std::vector<float> upper_dist_matrix;   //CalcLJPotential
    std::vector<float> pair_index_a;    //TermCalcLJPotential
    std::vector<float> pair_index_b;    //TermCalcLJPotential
    std::vector<float> pairwise_dist_matrix;    //TermCalcLJPotential
    std::vector<std::vector<float>> atom_pairwise_dist_matrix;  //TermCalcLJPotential
    std::vector<float> lj_energy;   // TermCalcLJPotential CalcLJPotential
    std::vector<std::vector<float>> atom_pairwise_dist_lj_energy_matrix;    //TermCalcLJPotential
    float sum_lj_energy;    //TermCalcLJPotential
    std::vector<float> potential_energy;

    //FORCE Var//
    std::vector<std::vector<float>> atom_pairwise_dist_lj_force_matrix;
    std::vector<float> lj_force;
    std::vector<float> forces_vector;
    std::vector<float> total_forces_vec_atoms;

private:
    float crack;

};




#endif
