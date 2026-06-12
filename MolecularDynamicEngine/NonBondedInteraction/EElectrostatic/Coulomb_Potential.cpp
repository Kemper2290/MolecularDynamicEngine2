

#include "Coulomb_Potential.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "../../PDB_Input/pdbreader.h"
#include "../../Utility/UtilityFunctions.h"


CoulombPotential::CoulombPotential()
{

}


float CoulombPotential::CalcCoulombPotential(std::vector<std::vector<float>> position_inp,const float& atom_charge, bool pbc, float boxlength)
{
    std::cout<< "Beggining of CalcCoulombPotential " << '\n';

    if (pbc == true)
    {
        std::cout << "PBC CalcCoulombPotential is turned on" << '\n';
    }
    else
    {
        std::cout << "PBC CalcCoulombPotential is turned off" << '\n';
    }

    std::vector<float> last_row = position_inp.back();

    int N_atoms = last_row.size()/3;
    float U_total = 0.0;

    for (int i =0 ; i < N_atoms;++i)
    {
        float xi = last_row[(i*3)];
        float yi = last_row[(i*3)+1];
        float zi = last_row[(i*3)+2];

        for (int j=i+1;j < N_atoms; ++j)
        {
            float xj = last_row[(j*3)];
            float yj = last_row[(j*3)+1];
            float zj = last_row[(j*3)+2];

            // displacement
            float dx = xi - xj;
            float dy = yi - yj;
            float dz = zi - zj;

            if (pbc == true)
            {
                dx -= boxlength * std::round(dx/boxlength);
                dy -= boxlength * std::round(dy/boxlength);
                dz -= boxlength * std::round(dz/boxlength);
            }

            else if (pbc == false)
            {

            }

            float r = std::sqrt(((dx*dx) + (dy*dy) + (dz*dz)));
            std::cout << "r value of coulomb: :" << r << "\n";
            float inv_r = 1.0f/r;

            if (r < 1e-12f) continue;

            const float ECfactor = 138.935458;
            float Ucp = ECfactor * atom_charge*atom_charge*inv_r;

            U_total += Ucp;

        }
    }

    coulomb_potential_energy.push_back(U_total);
    std::cout << "End of CalcCoulombPotential" << '\n';
    return U_total;

}

std::vector<float> CoulombPotential::CalcCoulombForce(std::vector<std::vector<float>> position_inp,const float& atomcharge,bool pbc,float boxlength)
{
    std::cout << "Beggining of CalcCoulombForce " << '\n';
    if (pbc == true)
    {
        std::cout << "PBC is turned on" << '\n';
    }
    else
    {
        std::cout << "PBC is turned off" << '\n';
    }

    std::vector<float> last_row = position_inp.back();
    int N_atoms = last_row.size()/3;
    std::vector<float> forces(last_row.size(),0);

    for (int i =0; i < N_atoms;++i)
    {
        float xi = last_row[(i*3)];
        float yi = last_row[(i*3)+1];
        float zi = last_row[(i*3)+2];

        for (int j= i + 1; j < N_atoms; ++j)
        {
            float xj = last_row[(j*3)];
            float yj = last_row[(j*3)+1];
            float zj = last_row[(j*3)+2];

            // displacement dx dy dz

            float dx = xi - xj;
            float dy = yi - yj;
            float dz = zi - zj;

            if (pbc == true)
            {
                dx -= boxlength * std::round(dx/boxlength);
                dy -= boxlength * std::round(dy/boxlength);
                dz -= boxlength * std::round(dz/boxlength);

            }
            else if (pbc == false)
            {
                //std::cout << "PBC force is off" << '\n';
            }


            float r2 = (dx*dx) + (dy*dy) + (dz*dz);
            float inv_r2 = 1.0f/r2;

            float force_scalar =138.935458 * inv_r2*atomcharge*atomcharge;
            std::cout << "Coulomb force_sclar: " << force_scalar << '\n';

            float fx = force_scalar * dx;
            float fy = force_scalar * dy;
            float fz = force_scalar * dz;

            // newtons 3rd law applied

            forces[(3*i)] += fx;
            forces[(3*i)+1] += fy;
            forces[(3*i)+2] += fz;

            forces[(3*j)] -= fx;
            forces[(3*j)+1] -= fy;
            forces[(3*j)+2] -= fz;
        }
    }

    std::cout <<'\n';
    //std::cout << "Forces vector fx1 fy1 fz1 fx2 " << '\n';
    //Print1DVec(forces);
    std::cout << " end of Calc Coulomb Force " << '\n';

    return forces;
}