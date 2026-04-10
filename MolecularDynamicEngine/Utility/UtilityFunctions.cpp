#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <variant>
#include <random>
#include <cmath>
#include "../matplotlibcpp.h"

#include "UtilityFunctions.h"

template <typename T>
void Print1DVec(T vec)
{
    for (auto& element: vec)
    {
        std::cout << element << " ";
    }
    std::cout << std::endl;
}

template void Print1DVec<std::vector<float>>(std::vector<float>);
template void Print1DVec<std::vector<int>>(std::vector<int>);
template void Print1DVec<std::vector<std::string>>(std::vector<std::string>);

////////////////////////
template <typename T>
void Print2DVec(T vec)
{
    for (auto& row: vec)
    {
        for (auto& element: row)
        {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template void Print2DVec<std::vector<std::vector<float>>>(std::vector<std::vector<float>>);

/////////////////////////////////////////////////
float CalcDist3d(float x1, float y1, float z1, float x2, float y2, float z2)
{   // dist = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)

    float dx = x2-x1;
    float dy = y2-y1;
    float dz = z2-z1;

    return std::sqrt((dx*dx) + (dy*dy) + (dz*dz));

}

float PBCCalcDist3d(float x1, float y1, float z1, float x2, float y2, float z2, float box_length)
{
    float dx = x1 - x2;
    dx -= box_length * std::round(dx/box_length);

    float dy = y1 - y2;
    dy -= box_length * std::round(dy/box_length);

    float dz = z1 - z2;
    dz -= box_length * std::round(dz/box_length);

    float r = std::sqrt((dx*dx) + (dy*dy) + (dz*dz));
    std::cout << "r distance: " << r << std::endl;
    return r;
}


std::vector<float> Gen_Velocity(std::vector<std::vector<float>> pdb2gmx_coord_input,
                                                float mass_input,
                                                float temp_k)
{
    std::cout << " Beginning of Gen_velocity" << std::endl;

    std::vector<float> last_row = pdb2gmx_coord_input.back();

    int N = (last_row.size()/3); // number of particles
    float T = 300.0; // temperature in kelvin
    float kB = 0.0083144621; // boltzmann constant in kJ/mol*k
    float mass = mass_input; // g/mol


    //standard deviation in nm/ps for random generation around mean
    float sigma = std::sqrt((kB * T) / mass);

    std::vector<float> vx(N), vy(N), vz(N);

    // generate generator
    std::mt19937 rng(42);
    std::normal_distribution<float> dist(0.0,sigma);

    // Generate Velocities
    for (int i = 0; i < N; ++i)
    {
        vx[i] = dist(rng);
        vy[i] = dist(rng);
        vz[i] = dist(rng);
    }

    //std::cout << "print vx " << std::endl;
    //Print1DVec(vx);

    //remove center-of-mass velocity
    float vx_com = 0.0, vy_com = 0.0, vz_com = 0.0;

    for (int i =0 ; i < N; ++i)
    {
        vx_com += vx[i];
        vy_com += vy[i];
        vz_com += vz[i];
    }

    vx_com /= N;
    vy_com /= N;
    vz_com /= N;

    //subtract com from vx to negate bulk momentum
    for (int i = 0; i < N; ++i)
    {
        vx[i] -=vx_com;
        vy[i] -=vy_com;
        vz[i] -=vz_com;
    }

    //compute current temperature for scaling since were taking a distribution that is randomly generate
    float kinetic = 0.0;
    for (int i = 0; i < N; ++i)
    {
        kinetic += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }

    float currentTemp = (mass * kinetic) / (3.0 * N * kB);

    // -- rescale velocities --
    float scale = std::sqrt(T / currentTemp);

    for (int i = 0; i < N; ++i)
    {
        vx[i] *= scale;
        vy[i] *= scale;
        vz[i] *= scale;
    }



    std::cout << " Vx,Vy,Vz after com removal and temp scale: " << std::endl;
    std::cout << "Vx: ";
    Print1DVec(vx);
    std::cout << "Vy: " ;
    Print1DVec(vy);
    std::cout << "Vz: " ;
    Print1DVec(vz);

    std::vector<float> gen_velocity_vecxyz;

    for (int i = 0; i < N; ++i)
    {
        gen_velocity_vecxyz.push_back(vx[i]);
        gen_velocity_vecxyz.push_back(vy[i]);
        gen_velocity_vecxyz.push_back(vz[i]);

    }

    Print1DVec(gen_velocity_vecxyz);
    std::cout << "End of gen_Velocity" << std::endl;
    return gen_velocity_vecxyz;
}


void Plot_Trajectory(std::vector<std::vector<float>> trajectory_input,float delay, float xlim,float ylim)
{
    // trajectory_input is a 2d matrix in format of {t=0{x1,y1,z1,x2,y2....},
    //                                               t=1{x1,y1,z1,x2,y2....}}

    std::vector<std::vector<float>> xcol;
    std::vector<std::vector<float>> ycol;
    std::vector<std::vector<float>> zcol;

    for (int i =0; i < trajectory_input.size(); ++i)
    {
        std::vector<float> xcol_vec;
        std::vector<float> ycol_vec;
        std::vector<float> zcol_vec;
        for (int j = 0; j < trajectory_input[0].size(); j+=3)
        {

            float xvalue = trajectory_input[i][j];
            float yvalue = trajectory_input[i][j+1];
            float zvalue = trajectory_input[i][j+2];

            xcol_vec.push_back(xvalue);
            ycol_vec.push_back(yvalue);
            zcol_vec.push_back(zvalue);
        }
        xcol.push_back(xcol_vec);
        ycol.push_back(ycol_vec);
        zcol.push_back(zcol_vec);
    }

    //std::cout << "Printing 2d vectors: xcol, ycol,zcol"<< std::endl;
    //Print2DVec(xcol);
    //std::cout << std::endl;
    //Print2DVec(ycol);
    //std::cout << std::endl;
    //Print2DVec(zcol);
    //std::cout << std::endl;


    namespace plt= matplotlibcpp;

    for (int frame =0; frame < trajectory_input.size();++frame)
    {

        //break;
        plt::clf();
        //std::string label = "Simulation Time: " + std::to_string(frame) + " Simulation Temp: " + std::to_string(temperature_list[frame]) + "F";
        std::string label = "Simulation Time: " + std::to_string(frame) + " Simulation Temp: 300K ";
        plt::named_plot(label,xcol[frame],ycol[frame],"bo");
        //plt::plot(xcol[0],ycol[0],"ro");
        //plt::plot(flower_xcol,flower_ycol,"go");
        plt::xlim(-xlim,ylim);
        plt::ylim(-ylim,ylim);
        plt::legend();
        plt::draw();
        plt::pause(delay);
    }
}

void Plot2D(std::vector<std::vector<float>> trajectory_input,float delay, float xlim,float ylim)
{
    std::vector<float> xcol;
    std::vector<float> ycol;

    for (int i =0; i < trajectory_input[0].size(); i+=3)
    {

        float x = trajectory_input[0][i];
        float y = trajectory_input[0][i+1];
        std::cout << "x: " << x << " y: " << y << std::endl;

        xcol.push_back(x);
        ycol.push_back(y);
    }

    std::cout << "xcol.size(): " << xcol.size() << std::endl;
    std::cout << "ycol.size(): " << ycol.size() << std::endl;

    namespace plt = matplotlibcpp;


    plt::plot(xcol,ycol,"ro");
    //plt::plot(xcol[0],ycol[0],"ro");
    //plt::plot(flower_xcol,flower_ycol,"go");
    plt::xlim(-xlim,ylim);
    plt::ylim(-ylim,ylim);
    //plt::legend();
    //plt::draw();
    plt::pause(delay);


}


// Apply PBC here did not work great because it used absolute position instead of relative position which could cause problems
/*
std::vector<float> ApplyPeriodicBoundaryConditions(std::vector<float> old_position_input,float box_length)
{
    std::cout << "Beginning of ApplyPeriodicBoundaryConditions: " << std::endl;

    std::vector<float> new_pbc_atom_position;

    for (int i = 0; i < old_position_input.size(); ++i)
    {
        std::cout << "old_position_input[" << i << "] = " << old_position_input[i] << std::endl;
        float floor = std::floor(old_position_input[i]/box_length);
        std::cout << "floor value: " << floor << std::endl;
        float lfloor = box_length * floor;
        std::cout << "lfloor value: " << lfloor << std::endl;
        float position = old_position_input[i] - lfloor;
        std::cout << "position after boxlength shit : " << position << std::endl;
        new_pbc_atom_position.push_back(position);
    }

    std::cout << " End of ApplyPeriodicBoundaryConditions: " << std::endl;
    return new_pbc_atom_position;
}
*/


std::vector<std::vector<float>> CreateAtomCluster(float spacing, float gridsize)
{
    std::vector<std::vector<float>> coordinates;

    for (int i = 0; i < gridsize; ++i)
    {
        for (int j = 0; j< gridsize; ++j)
        {
            for (int k =0; k < gridsize; ++k)
            {
                float x_coord = i * spacing;
                float y_coord = j * spacing;
                float z_coord = k * spacing;

                coordinates.push_back({x_coord, y_coord, z_coord});
            }
        }
    }

    Print2DVec(coordinates);

    return coordinates;
}












