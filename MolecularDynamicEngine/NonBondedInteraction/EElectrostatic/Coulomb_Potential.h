#ifndef COULOMB_POTENTIAL_H
#define COULOMB_POTENTIAL_H

#include <string>
#include <vector>
#include <cmath>
#include <iostream>

class CoulombPotential
{
public:
    CoulombPotential();

    float CalcCoulombPotential(
        std::vector<std::vector<float>> position_inp,
        const float& atom_charge,
        bool pbc,
        float boxlength
        );

    std::vector<float> CalcCoulombForce(
        std::vector<std::vector<float>> position_inp,
        const float& atom_charge,
        bool pbc,
        float boxlength
        );


    //Variables
    std::vector<float> coulomb_potential_energy;


private:

};



#endif





