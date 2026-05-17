
#include "Lennard_Jones_Potential.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "../../PDB_Input/pdbreader.h"
#include "../../Utility/UtilityFunctions.h"
#include "../../Movement/VelocityPositionAlgo.h"



LJPotential::LJPotential()
{

}



void LJPotential::ReadLJParameters(const std::string& filename)
{
    std::cout << "Reading the Lennard-Jones Parameters.txt file "<<std::endl;
    std::cout << std::endl;

    std::ifstream file(filename);// opens the file filename. file is the variable name
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;// if file doesnt open this error occurs
        return;
    }
    std::string line; // variable to store each line from file
    while (std::getline(file,line))
    {
        if (line.substr(0,4)=="PAIR")
        {
            std::string atom_pair = line.substr(8,5);
            lj_atompair.push_back(atom_pair);
            float sigma = std::stof(line.substr(21,4));
            sigma = sigma * std::pow(10,-1);
            lj_sigma.push_back(sigma);
            float epsilon = std::stof(line.substr(37,4));
            lj_epsilon.push_back(epsilon);

            std::cout << "Atom pair: " << atom_pair
                      << ", Sigma Value: " << sigma
                      << ", Epsilon Value: " << epsilon << std::endl;
            std::cout << std::endl;
        }
    }
}


float LJPotential::CalcLJPotential(
                        std::vector<std::vector<float>> pdb2gmx_coord_input,
                        std::vector<std::vector<float>>dist_matrix_input,
                        std::vector<std::string> pdb2gmx_atomlist_input)
{
    // input the pdb2gmx gmx input which should be a 2d matrix of form {{x1,y1,z1,x2,y2,z2,..}}
    //      and use the end row to calculate the potential energy
    //      a seperate function will calculate the force generated from it

    std::cout << "Beginning of CalcLJPotential Function: " << std::endl;
    std::cout << std::endl;


    std::cout << "Printing dist_matrix_input " << std::endl;
    Print2DVec(dist_matrix_input);

    // Loop takes the dist_matrix_input and grabs only the upper triangle matrix and stores in a 1d vector
    for (int i=0;i < dist_matrix_input.size();++i)
    {
        for (int j=0; j < dist_matrix_input[0].size();++j)
        {
            if (j>=i)
            {
                upper_dist_matrix.push_back(dist_matrix_input[i][j]);
            }
        }
    }
    std::cout << std::endl;
    std::cout <<"Printing upper_dist_matrix: " << std::endl;
    Print1DVec(upper_dist_matrix);
    std::cout << std::endl;

    float fourep= 4* lj_epsilon[0];

    // V(r) = 4e * ((sigma/r)^12 + (sigma/r)^6)
    for (int i=0; i < upper_dist_matrix.size(); ++i)
    {
        if (upper_dist_matrix[i] != 0.0f)
        {
            std::cout << "Upper_dist_matrix[i] != 0" << std::endl;
            float sig12 = pow((lj_sigma[0]/upper_dist_matrix[i]),12);
            float sig6 = pow((lj_sigma[0]/upper_dist_matrix[i]),6);
            std::cout << "sig12: " << sig12 << " sig6: " << sig6 << std::endl;

            float energy = fourep * (sig12 - sig6);
            lj_energy.push_back(energy);

        }
        else if (upper_dist_matrix[0] == 0.0f)
        {
            std::cout << "Upper_dist_matrix[i] = 0" << std::endl;
            float sig12 = 0;
            float sig6 = 0;
            std::cout << "sig12: " << sig12 << " sig6: " << sig6 << std::endl;

            float energy = 0;
            lj_energy.push_back(energy);

        }
    }
    std::cout << std::endl;
    std::cout << "lj_energy: ";
    Print1DVec(lj_energy);

    sum_lj_energy=0;
    for (int i=0;i < lj_energy.size(); ++i)
    {
        sum_lj_energy += lj_energy[i];
    }
    std::cout << "Sum_lj_energy: " << sum_lj_energy << std::endl;


    std::cout << "END CalcLJPotential Function:" << std::endl;
    return sum_lj_energy;

}

// i believe i need to calculate the force for each atom so the above method doesn't work correctly i believe
//      since i need to eventually get the velcoity of each atom also

std::vector<std::vector<float>> LJPotential::PairwiseCalcLJPotential(
                        std::vector<std::vector<float>> pdb2gmx_coord_input,
                        std::vector<std::string> pdb2gmx_atomlist_input)
{
    std::cout << "Beginning of CalcLJPotential Function: " << std::endl;
    std::cout << std::endl;


    pair_index_a.clear();
    pair_index_b.clear();
    pairwise_dist_matrix.clear();
    atom_pairwise_dist_matrix.clear();
    atom_pairwise_dist_lj_energy_matrix.clear();



    std::vector<float> last_row = pdb2gmx_coord_input.back();
    std::cout << "last row size: " << last_row.size() << std::endl;
    std::cout << std::endl;
    int a_increase = 1;

    // ////////////// BEGIN OF ATOM INDEX PAIR WISE DISTANCE MATRIX SETUP ////////////////////////////////////////
    for (int i =0; i < last_row.size(); i+= 3) // part of creating distance matrix that shows pairwise distance between atoms
    {
        pair_index_a.insert(pair_index_a.end(),(last_row.size()/3),a_increase);
        a_increase = a_increase + 1;
    }

    for (int i=0;i < last_row.size(); i +=3)
    {
        int b_increase = 1;
        for (int j = 0; j < (last_row.size()/3); ++j)
        {
            pair_index_b.push_back(b_increase);
            b_increase = b_increase + 1;
        }
    }
    // ////////// END OF ATOM INDEX PAIRWISE DISTANCE MATRIX SETUP ///////////////////////////

    // Beginning of LJ potential energy calculation

    float fourep= 4* lj_epsilon[0];
    for (int i =0; i < last_row.size(); i+=3)
    {
        for (int j = 0; j < last_row.size(); j+=3)
        {
            float result = CalcDist3d(
                last_row[i],last_row[i+1],last_row[i+2],
                last_row[j],last_row[j+1],last_row[j+2]);

            //std::cout << "These are the result: " << result << std::endl;
            pairwise_dist_matrix.push_back(result);

            if (result == 0.0f)
            {
                float sig12 = 0;
                float sig6 = 0;
                float energy = 0;
                //std::cout << "sig12: " << sig12 << " sig6: " << sig6 << "energy: " << energy << std::endl;
                lj_energy.push_back(0); // could have issues when running the simulation if variables dont reset
            }
            else if (result != 0.0f)
            {
                // change to useing sig*sig method to speed up calculation
                float sig12 = pow((lj_sigma[0]/result),12);
                float sig6 = pow((lj_sigma[0]/result),6);
                float energy = fourep * (sig12-sig6);
                //std::cout << "sig12: " << sig12 << " sig6: " << sig6 << " energy: " << energy << std::endl;
                lj_energy.push_back(energy);
            }
        }
    }


    // I fucked up the loop and the matrix isnt the right size
    for (int i = 0; i < pair_index_a.size(); ++i)
    {
        atom_pairwise_dist_matrix.push_back({pair_index_a[i],
                                                   pair_index_b[i],
                                                   pairwise_dist_matrix[i]});
    }

    std::cout << std::endl;
    //std::cout << "Pair_index_a vector: " << std::endl;
    //Print1DVec(pair_index_a);
    //std::cout << std::endl;
    //std::cout << "Pair_index_b vector: " << std::endl;
    //Print1DVec(pair_index_b);
    //std::cout << "atom_pairwise_dist_matrix: Atom Index starts at 1 " << std::endl;
    //Print2DVec(atom_pairwise_dist_matrix);

    std::cout << std::endl;
    std::cout << "lj_energy: ";
    Print1DVec(lj_energy);

    for (int i = 0; i < pair_index_a.size();++i)
    {
        atom_pairwise_dist_lj_energy_matrix.push_back({pair_index_a[i],
                                                        pair_index_b[i],
                                                        pairwise_dist_matrix[i],
                                                        lj_energy[i]});
    }
    std::cout << "atom_pairwise_lj_energy_matrix: Atom Index starts at 1: " << std::endl;
    std::cout << "Index Index Dist Energy(kJ/mol) " << std::endl;
    Print2DVec(atom_pairwise_dist_lj_energy_matrix);

    sum_lj_energy = 0;
    for (int i = 0; i < lj_energy.size(); ++i)
    {
        sum_lj_energy = sum_lj_energy + lj_energy[i];
    }
    sum_lj_energy = sum_lj_energy / 2;
    std::cout << "sum_lj_energy: " << sum_lj_energy << std::endl;


    a_increase = 1; // hopefully reset var if function is called again


    std::cout << std::endl;
    std::cout << "END TermCalcLJPotential Function:" << std::endl;
    std::cout << std::endl;

    //std::vector<std::vector<float>> trash;
    return atom_pairwise_dist_lj_energy_matrix;
}

float LJPotential::PBCCalcLJPotential(std::vector<std::vector<float>> position_inp,
                                        float sigma,
                                        float epsilon,
                                        bool pbc,
                                        float box_length)
{
    std::cout << "Beginning of PBCCalcLJPotential " << std::endl;

    if (pbc == true)
    {
        std::cout << "PBC LJ Potential is turned on " << std::endl;
    }
    else
    {
        std::cout << "PBC LJ Potential is turned off " << std::endl;
    }


    std::vector<float> last_row = position_inp.back();

    int N_atoms = last_row.size()/3;
    float U_total = 0.0;

    float sigma6 = pow(sigma,6);
    float sigma12 = sigma6 * sigma6;

    for (int i = 0; i < N_atoms; ++i)
    {
        float xi = last_row[3*i];
        float yi = last_row[(3*i) + 1];
        float zi = last_row[(3*i) + 2];

        for (int j = i + 1; j < N_atoms; ++j)
        {
            float xj = last_row[3*j];
            float yj = last_row[3*j + 1];
            float zj = last_row[3*j + 2];

            // --- displacement ---
            float dx = xi - xj;
            float dy = yi - yj;
            float dz = zi - zj;

            if (pbc == true)
            {
                dx -= box_length * std::round(dx/box_length);
                dy -= box_length * std::round(dy/box_length);
                dz -= box_length * std::round(dz/box_length);
            }

            else if (pbc == false)
            {

            }

            float r2 = ((dx*dx) + (dy*dy) + (dz*dz));

            // avoids division by zero and so distance to self isnt counted
            if (r2 < 1e-12f) continue;

            float inv_r2 = 1.0f/r2;
            float inv_r6 = inv_r2*inv_r2*inv_r2;
            float inv_r12 = inv_r6*inv_r6;

            // Lennard-Jones potential

            float Uij = 4 * epsilon * ((sigma12 * inv_r12) - (sigma6*inv_r6));
            U_total += Uij;

        }
    }

    potential_energy.push_back(U_total);
    std::cout << " End of PBCCalcLJPotential " << std::endl;
    return U_total;
}


std::vector<float> LJPotential::CalcForceLJPotential(
                        std::vector<std::vector<float>> pdb2gmx_coord_input,
                        std::vector<std::string> pdb2gmx_atomlist_input,
                        bool pbc,
                        float box_length_input)
{
    // Using equation F= ((48*e)/sigma)*(((sigma/r)^13) -(.5*(sigma/r)^7))
    std::cout << "Beginning of CalcForceLJPotential Function: " << std::endl;
    std::cout << std::endl;



    pair_index_a.clear();
    pair_index_b.clear();
    pairwise_dist_matrix.clear();
    atom_pairwise_dist_matrix.clear();
    total_forces_vec_atoms.clear();
    lj_force.clear();
    forces_vector.clear();
    atom_pairwise_dist_lj_force_matrix.clear();

    std::vector<float> last_row = pdb2gmx_coord_input.back();
    std::cout << "last row size: " << last_row.size() << std::endl;
    std::cout << std::endl;

    /* // PBC not working great here
    if (pbc==true)
    {
        std::cout << "PBC in CalcForceLJPotential is enabled " << std::endl;
        last_row = ApplyPeriodicBoundaryConditions(last_row,box_length_input);

    }
    else
    {
        std::cout << "PBC in CalcForceLJPotential is not enables " << std::endl;
    }
    */

    // ////////////// BEGIN OF ATOM INDEX PAIR WISE DISTANCE MATRIX SETUP ////////////////////////////////////////
    int a_increase = 1;
    for (int i =0; i < last_row.size();i+=3)
    {
        pair_index_a.insert(pair_index_a.end(),(last_row.size()/3),a_increase);
        a_increase = a_increase + 1;
    }

    for (int i =0; i < last_row.size(); i +=3)
    {
        int b_increase = 1;
        for (int j = 0; j< (last_row.size()/3); ++j)
        {
            pair_index_b.push_back(b_increase);
            b_increase = b_increase + 1;
        }
    }

    // ////////// END OF ATOM INDEX PAIRWISE DISTANCE MATRIX SETUP ///////////////////////////

    // Beginning of LJ Net Force

    float four8epsig = ((48*lj_epsilon[0])/lj_sigma[0]);
    //std::cout << "four8epsig: " << four8epsig << std::endl;
    for (int i = 0; i < last_row.size();i+=3)
    {
        for (int j = 0; j<last_row.size();j+=3)
        {
            float result = CalcDist3d(
                    last_row[i], last_row[i+1], last_row[i+2],
                    last_row[j], last_row[j+1], last_row[j+2]);


            float inv_result = 1 / result;
            //std::cout << "These are the result: " << result << std::endl;
            pairwise_dist_matrix.push_back(result);

            if (result == 0.0f)
            {
                float sig13 = 0;
                float sig7 = 0;
                float force = 0;
                //std::cout << "sig13: " << sig13 << " sig7: " << sig7 << "Total force: " << force << std::endl;

                lj_force.push_back(0);

                float dx = last_row[j] - last_row[i];
                float dy = last_row[j+1] - last_row[i+1];
                float dz = last_row[j+2] - last_row[i+2];

                float force_vecx = 0;
                float force_vecy = 0;
                float force_vecz = 0;

                forces_vector.push_back(force_vecx);
                forces_vector.push_back(force_vecy);
                forces_vector.push_back(force_vecz);

            }
            else if (result != 0.0f)
            {
                float sig13 = pow((lj_sigma[0]/result),13);
                float sig7 = pow((lj_sigma[0]/result),7);
                float force = four8epsig*((sig13)-(.5*(sig7))); // this is magnitude
                //std::cout << "sig13: " << sig13 << " sig7: " << sig7 << "Total force: " << force << std::endl;
                lj_force.push_back(force);
                float dx = last_row[j] - last_row[i];
                float dy = last_row[j+1] - last_row[i+1];
                float dz = last_row[j+2] - last_row[i+2];

                float force_vecx = force * dx * inv_result;
                float force_vecy = force * dy * inv_result;
                float force_vecz = force * dz * inv_result;

                forces_vector.push_back(force_vecx);
                forces_vector.push_back(force_vecy);
                forces_vector.push_back(force_vecz);

            }
        }
    }

    // creates this matrix that contains atom indexs, distance, ||force||
    for (int i = 0; i < pair_index_a.size(); ++i)
    {
        atom_pairwise_dist_matrix.push_back({pair_index_a[i],
                                                   pair_index_b[i],
                                                   pairwise_dist_matrix[i]});
    }

    for (int i = 0; i < pair_index_a.size();++i)
    {
        atom_pairwise_dist_lj_force_matrix.push_back(
                                          {pair_index_a[i],
                                                pair_index_b[i],
                                                pairwise_dist_matrix[i],
                                                lj_force[i]});
    }

    // perform vector addition for forces_vector to get total force vector on each atom

    total_forces_vec_atoms.resize(last_row.size(),0);
    int a = 0;
    std::cout <<" int a = 0: " << a << std::endl;
    for (int i = 0 ; i < forces_vector.size(); i+=forces_vector.size()/3) //9// i+= forces_vector.size()/3 is same size as last_row
    {
        for (int j =i; j<i+(forces_vector.size()/3);j+=3 )
        {
            //std::cout << "j=1 is: " << j << std::endl;
            //std::cout << "forces_vector[" << j << "]: " << forces_vector[j] << std::endl;
            //std::cout << "forces_vector[" << j+1 << "]: " << forces_vector[j+1] << std::endl;
            //std::cout << "forces_vector[" << j+2 << "]: " << forces_vector[j+2] << std::endl;

            total_forces_vec_atoms[a+0] += forces_vector[j];
            total_forces_vec_atoms[a+1] += forces_vector[j+1];
            total_forces_vec_atoms[a+2] += forces_vector[j+2];
        }
        a+=3;
        //std::cout << "a+=3: " << a << std::endl;
    }


    ///*
    //std::cout << "Pair_index_a vector: " << std::endl;
    //Print1DVec(pair_index_a);
    //std::cout << std::endl;
    //std::cout << "Pair_index_b vector: " << std::endl;
    //Print1DVec(pair_index_b);
    //std::cout << "atom_pairwise_dist_matrix: Atom Index starts at 1 " << std::endl;
    //Print2DVec(atom_pairwise_dist_matrix);
//
    //std::cout << std::endl;
    //std::cout << "lj_force: ";
    //Print1DVec(lj_force);
    //std::cout << std::endl;
//
//
    //std::cout << "atom_pairwise_dist_lj_force_matrix: Atom Index starts at 1:" << std::endl;
    //std::cout << "Index Index Dist(A) ||Force|| (kJ/(mol*A) " << std::endl;
    //Print2DVec(atom_pairwise_dist_lj_force_matrix);
    //*/

    std::cout << "Forces_vector: fx1,fy1,fz1,fx2,fy2,fz2,fx3... " << std::endl;
    Print1DVec(forces_vector);
    std::cout << "forces_vector.size(): " << forces_vector.size() << std::endl;

    std::cout << std::endl;
    std::cout << "Total_forces_vec_atoms: total force vectors for each atom from vector addition "<< std::endl;
    Print1DVec(total_forces_vec_atoms);


    a_increase = 1;
    std::cout << std::endl;
    std::cout << "END CalcForceLJPOTENTIAL Function: " << std::endl;
    std::cout << std::endl;

    crack= 987;
    return total_forces_vec_atoms; //1d vector

}




std::vector<float> LJPotential::PBCForceLJPotential(
                        std::vector<std::vector<float>> position_inp,
                        float sigma, float epsilon,
                        bool pbc,float box_length)
{
    // using form F = 24epsilon * (2sigma12/r14 - sigma6/r8)
    std::cout << " Begin PBCForceLJPotential " << std::endl;
    if (pbc == true)
    {
        std::cout << "PBC Force is turned on " << std::endl;
    }
    else
    {
        std::cout << "PBC Force is turned off " << std::endl;
    }

    std::vector<float> last_row = position_inp.back();

    int N_atoms = last_row.size()/3;
    std::vector<float> forces(last_row.size(),0);
    float sigma6 = std::pow(sigma,6);
    float sigma12 = sigma6 * sigma6;

    for (int i =0; i < N_atoms; ++i)
    {
        //std::cout << "force i : " << i << std::endl;
        float xi = last_row[3*i];
        float yi = last_row[(3*i)+1];
        float zi = last_row[(3*i)+2];

        for ( int j = i +1; j < N_atoms;++j)
        {
            //std::cout << "force j : " << j << std::endl;
            float xj = last_row[3*j];
            float yj = last_row[(3*j)+1];
            float zj = last_row[(3*j)+2];

            // displacement
            float dx = xi - xj;
            float dy = yi - yj;
            float dz = zi - zj;

            if (pbc==true)
            {
                //std::cout << "PBC Force is turned on " << std::endl;
                // dx -= is the same as dx = dx - stuff but only evaluates once
                dx -= box_length * std::round(dx/box_length);   // round is used here instead of floor on purpose
                dy -= box_length * std::round(dy/box_length);
                dz -= box_length * std::round(dz/box_length);

            }
            else if (pbc==false)
            {
                //std::cout << "PBC Force is turned off ";
            }

            float r2 = dx*dx + dy*dy + dz*dz;
            //std::cout << "r2 is: " << r2 << std::endl;

            if (r2 < 1e-12) continue; // continues to next loop since r = 0 which means distance is so far its 0 because were only measuring unique pairs

            float inv_r2 = 1.0f / r2;
            float inv_r6 = inv_r2 * inv_r2 * inv_r2;
            float inv_r12 = inv_r6 * inv_r6;



            float force_scalar = 24*epsilon * ((2*sigma12*inv_r12) - (sigma6*inv_r6)) * inv_r2;

            float fx = force_scalar * dx;
            float fy = force_scalar * dy;
            float fz = force_scalar * dz;

            // Newtons 3rd law applied

            forces[3*i]       += fx;
            forces[(3*i) + 1] += fy;
            forces[(3*i) + 2] += fz;

            forces[3*j]       -= fx;
            forces[(3*j) + 1] -= fy;
            forces[(3*j) + 2] -= fz;

        }
    }
    std::cout << std::endl;
    std::cout << "Forces vectors: fx1 fy1 fz1 fx2 " << std::endl;
    Print1DVec(forces);
    std::cout << " End of PBCForceLJPotential " << std::endl;

    return forces;
}





