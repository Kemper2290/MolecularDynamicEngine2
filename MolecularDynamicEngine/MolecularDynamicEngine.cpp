#include <vector>
#include "PDB_Input/myheader.h" // Include your user-defined header
#include <iostream>
#include <algorithm>

#include "PDB_Input/pdbreader.h"
#include <Python.h>
#include "matplotlibcpp.h"
#include "NonBondedInteraction/EVanDerWaals/Lennard_Jones_Potential.h"
#include "Utility/UtilityFunctions.h"
#include "Movement/VelocityPositionAlgo.h"

namespace plt= matplotlibcpp;

float timestep = 0.002; // unit is ps
float simtime = 1000; // units in ps
float argon_mass = 39.948; //g/mol
float box_size = 3;

int main() {
    // Use the function and class declared in myheader.h and defined in myheader.cpp

    std::cout << "Hello World" << std::endl;
    myFunction();
    MyClass obj;
    obj.doSomething();

    MyPdbInfo mypdb;
    mypdb.ReadPDB("di_argon.pdb");
    //mypdb.ReadPDB("tri_argon.pdb");
    //mypdb.ReadPDB("10_argon.pdb");
    //mypdb.ReadPDB("multi_argon.pdb");


    //mypdb.MeasureBondDist(mypdb.pdbmatrix);
    std::vector<std::vector<float>> pdb2gmx_coord = mypdb.Pdb2Gmx(mypdb.pdbmatrix);
    //std::vector<std::vector<float>> atoms = CreateAtomCluster(4,4);
    //std::vector<std::vector<float>> pdb2gmx_coord = mypdb.Pdb2Gmx(atoms);
    //std::vector<std::string> pdb2gmx_atomlist(16,"AR");
    std::vector<std::string> pdb2gmx_atomlist = mypdb.pdbatomname;


    // Lennard-Jones section
    LJPotential lj_object;
    lj_object.ReadLJParameters("NonbondedInteraction/EVanDerWaals/LJParameters.txt");

    std::vector<float> generated_init_velocity = Gen_Velocity(pdb2gmx_coord,argon_mass,300);
    std::vector<float> velocity = generated_init_velocity;
    std::vector<std::vector<float>> trajectory = pdb2gmx_coord;

    lj_object.PairwiseCalcLJPotential(trajectory,pdb2gmx_atomlist);
    float sigma_ar = lj_object.lj_sigma[0];
    float epsilon_ar = lj_object.lj_epsilon[0];
    float lj_energy = lj_object.PBCCalcLJPotential(trajectory,sigma_ar,epsilon_ar,true,box_size);
    std::cout << "LJ comparison energy : " << lj_energy << std::endl;

    std::vector<float> force = lj_object.PBCForceLJPotential(trajectory,sigma_ar,epsilon_ar,true,box_size);
    std::vector<float> acceleration = Movement::CalcAcceleration(argon_mass,force);


    for (int i = 0; i < simtime; ++i)
    {
        std::cout << "STARTING SIM LOOP: " << i << std::endl;

        std::vector<float> newposition = Movement::PBCCalcVerletPosition(trajectory,velocity,
                                                                        acceleration,timestep,box_size,true);
        float lj_energy = lj_object.PBCCalcLJPotential(trajectory,sigma_ar,epsilon_ar,true,box_size);
        std::cout <<" lj_energy: " << lj_energy << std::endl;
        std::vector<float> newforce = lj_object.PBCForceLJPotential(trajectory,sigma_ar,epsilon_ar,true,box_size);
        std::vector<float> newacceleration = Movement::CalcAcceleration(argon_mass,newforce);
        std::vector<float> newvelocity = Movement::CalcVerletUpdateVelocity(velocity,acceleration,
                                                            newacceleration,timestep, argon_mass);
        //std::cout << "NEWVELOCITY.SIZE(): " << newvelocity.size() << std::endl;
        trajectory.push_back(newposition);

        //std::cout << "Printing Trajectory: " << std::endl;
        //Print2DVec(trajectory);


        //force = newforce;
        acceleration = newacceleration;
        velocity = newvelocity;

        std::cout << "END OF SIM LOOP: " << i << std::endl;

    }

    std::cout << "Final Position Result From Trajectory Matrix: " << std::endl;
    Print2DVec(trajectory);

    std::vector<float> trajectory_time;
    std::vector<float> atom_dist;

    int total_time = trajectory.size();
    for (int i =0;i < trajectory.size();++i)
    {
        trajectory_time.push_back(i);
        float dist = PBCCalcDist3d(trajectory[i][0],trajectory[i][1],trajectory[i][2],
                                    trajectory[i][3],trajectory[i][4],trajectory[i][5],box_size);
        atom_dist.push_back(dist);
    }





    // --- PLOTTING SECTION ---
    std::vector<float> PE_output = lj_object.potential_energy;
    std::cout << "PE_output.size(): " << PE_output.size() << std::endl;
    int PE_size = PE_output.size();
    std::vector<int> PE_index(PE_size);
    std::iota(PE_index.begin(),PE_index.end(),0);
    auto maxPE_output = std::max_element(PE_output.begin(),PE_output.end());
    float maxPE = *maxPE_output;
    auto minPE_output = std::min_element(PE_output.begin(),PE_output.end());
    float minPE = *minPE_output;
    std::cout << "maxPE: " << maxPE << " minPE: " << minPE << std::endl;




    std::vector<float> KE_output = Movement::kinetic_energy;
    std::cout << "KE_output.size(): " << KE_output.size() << std::endl;
    int KE_size = KE_output.size();
    std::vector<int> KE_index(KE_size); // creates vector of size KE_output.size() filled with 0's
    std::iota(KE_index.begin(),KE_index.end(),0);
    auto maxKE_output = std::max_element(KE_output.begin(),KE_output.end());
    float maxKE = *maxKE_output;
    auto minKE_output = std::min_element(KE_output.begin(),KE_output.end());
    float minKE = *minKE_output;
    std::cout << "maxKE: " << maxKE << " minKE: " << minKE << std::endl;



    std::vector<float> TE_output;
    for (int i = 0; i < KE_size; ++i)
    {
        float total_energy = KE_output[i] + PE_output[i];
        TE_output.push_back(total_energy);
    }
    int TE_size = TE_output.size();
    std::vector<int> TE_index(TE_size);
    std::iota(TE_index.begin(),TE_index.end(),0);



    plt::plot(PE_index,PE_output);
    plt::title("Potential Energy");
    plt::xlim(0,PE_size + 50);
    //plt::ylim(minPE-(minPE*.1f),(maxPE+(maxPE*.1f)));
    plt::pause(11);

    plt::plot(KE_index,KE_output);
    plt::title("Kinetic Energy");
    plt::xlim(0,(KE_size + 50));
    //plt::ylim(0.0f,maxKE);
    plt::pause(11);

    plt::plot(TE_index,TE_output);
    plt::title("Total Energy (kJ/mol)");
    plt::xlim(0,TE_size + 50);
    plt::pause(11);


    plt::plot(trajectory_time,atom_dist);
    plt::title("Atom Distance (nm) ");
    plt::xlim(0,total_time);
    plt::pause(11);





    Plot_Trajectory(trajectory,.0005,6,6);

    float test = mypdb.pdbx[0];
    std::cout << "THIS IS VAR TEST: " << test << std::endl;

    return 0;
}









/*
for (int frame =0; frame < data.size();++frame)
    {
        //break;
        plt::clf();
        std::string label = "Simulation Time: " + std::to_string(frame) + " Simulation Temp: " + std::to_string(temperature_list[frame]) + "F";
        plt::named_plot(label,xcol[frame],ycol[frame],"bo");
        plt::plot(xcol[0],ycol[0],"ro");
        plt::plot(flower_xcol,flower_ycol,"go");
        plt::xlim(-300,300);
        plt::ylim(-300,300);
        plt::legend();
        plt::draw();
        plt::pause(0.05);
    }
*/