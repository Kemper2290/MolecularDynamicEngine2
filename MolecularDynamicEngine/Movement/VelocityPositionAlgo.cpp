#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // Required for fixed-width parsing (optional, can also use substr)
#include "../Utility/UtilityFunctions.h"
#include "../Movement/VelocityPositionAlgo.h"


// three main methods for calculating velocity and position
// 1) Euler Method: v(t+Δt)=v(t)+a(t)Δt and r(t+Δt)=r(t)+v(t)Δt
// 2) Velocity Verlet Method: step1 update position r(t+Δt)=r(t)+v(t)Δt+(1/2)a(t)Δt2. Step2 compute new force-> acceleration
//              Step 3. Update Velocity v(t+Δt)=v(t)+(1/2)(a(t)+a(t+Δt))Δt
// 3) Leapfrog method: v(t+Δt/2)


//std::vector<float> Movement::force;
std::vector<float> Movement::atom_accel;
std::vector<float> avg_atom_accel;
std::vector<float> Movement::new_atom_position;
std::vector<float> Movement::new_atom_velocity;
//std::vector<float> Movement::new_pbc_atom_position;
std::vector<float> Movement::kinetic_energy;
std::vector<float> Movement::momentum_magnitude;

std::vector<float> Movement::CalcAcceleration(float mass_input,
                                              std::vector<float> atom_force_input
)
{
    std::cout << " Beggining of Movement::CalcAcceleration " << std::endl;

    atom_accel.clear();

    for (int i = 0; i < atom_force_input.size(); ++i)
    {
        float ai = atom_force_input[i] / mass_input;
        atom_accel.push_back(ai);
    }

    std::cout << "atom_accel vectors: " << std::endl;
    Print1DVec(atom_accel);


    std::cout << "End of Movement::CalcACCeleration " << std::endl;
    std::cout << std::endl;
    return atom_accel;
}

std::vector<float> Movement::CalcVerletPosition(std::vector<std::vector<float>> coord_input,
                                                std::vector<float> atom_accel_input,
                                                std::vector<float> velocity_input,
                                                float dt,
                                                bool pbc,
                                                float box_length_input)
{
    // r(t+Δt)=r(t)+v(t)Δt+.5a(t)Δt2 utilizing verlet velocity position update
    // use xnew = =xold(t)+vx(t)Δt+.5ax(t)Δt2

    std::cout << std::endl;
    std::cout << " Beginning of CalcVerletPosition: " << std::endl;

    new_atom_position.clear();

    std::vector<float> last_row = coord_input.back();
    std::cout << "last_row coord: ";
    Print1DVec(last_row);

    if (pbc == true)
    {
        std::cout << "PBC in CalcVerletPosition is enabled " << std::endl;
        last_row = coord_input.back();
    }
    else
    {
        std::cout << "PBC in CalcVerletPosition is not enabled " << std::endl;
    }


    for (int i = 0; i < last_row.size(); i += 3)
    {
        float xold = last_row[i];
        float yold = last_row[i + 1];
        float zold = last_row[i + 2];

        float vxold = velocity_input[i];
        float vyold = velocity_input[i + 1];
        float vzold = velocity_input[i + 2];

        float vxdt = vxold * dt;
        float vydt = vyold * dt;
        float vzdt = vzold * dt;

        float axold = atom_accel_input[i];
        float ayold = atom_accel_input[i + 1];
        float azold = atom_accel_input[i + 2];

        float halfxaccel_dt = (1 / 2) * axold * dt * dt;
        float halfyaccel_dt = (1 / 2) * ayold * dt * dt;
        float halfzaccel_dt = (1 / 2) * azold * dt * dt;

        float xnew = xold + vxdt + halfxaccel_dt;
        float ynew = yold + vydt + halfyaccel_dt;
        float znew = zold + vzdt + halfzaccel_dt;

        new_atom_position.push_back(xnew);
        new_atom_position.push_back(ynew);
        new_atom_position.push_back(znew);
    }

    std::cout << "atom_new_position: " << std::endl;
    Print1DVec(new_atom_position);

    std::cout << std::endl;
    std::cout << "End of CalcVerletPosition: " << std::endl;
    std::cout << std::endl;
    return new_atom_position;
}

std::vector<float> Movement::PBCCalcVerletPosition(std::vector<std::vector<float>> position_inp,
                                                   std::vector<float> velocity_inp,
                                                   std::vector<float> acceleration_inp,
                                                   float dt,
                                                   float box_length,
                                                   bool pbc)
{
    std::cout << "Beginning of PBCCalcVerletPosition: " << std::endl;
    if (pbc == true)
    {
        std::cout << "PBC Verlet Position is turned on " << std::endl;
    }
    else
    {
        std::cout << "PBC Verlet Position is turned off " << std::endl;
    }


    std::vector<float> last_row = position_inp.back();
    std::vector<float> r_new(last_row.size());

    int N_atoms = last_row.size() / 3;

    for (int i = 0; i < N_atoms; ++i)
    {
        for (int d = 0; d < 3; ++d)
        {
            int idx = (3 * i) + d;

            r_new[idx] = last_row[idx] + velocity_inp[idx] * dt + 0.5 * acceleration_inp[idx] * dt * dt;
            if (pbc == true)
            {
                //std::cout << "PBC Verlet Position is turned on " << std::endl;
                r_new[idx] -= box_length * std::floor(r_new[idx] / box_length);
            }
            else if (pbc == false)
            {
                //std::cout << "PBC Verlet Position is turned off " << std::endl;
            }
        }
    }
    std::cout << "r_new position vector: " << std::endl;
    Print1DVec(r_new);

    std::cout << "End of PBCCalcVerletPosition: " << std::endl;

    return r_new;
}


std::vector<float> Movement::CalcVerletUpdateVelocity(std::vector<float> old_velocity_input,
                                                      std::vector<float> old_atom_accel_input,
                                                      std::vector<float> new_atom_accel_input,
                                                      float dt,
                                                      float mass_input)
{
    std::cout << "Beginning of CalcVerletUpdateVelocity: " << std::endl;
    //v(t+Δt) = v(t) + .5(a(t) + a(t+Δt))Δt

    new_atom_velocity.clear();


    float KE_total = 0.0f;
    float Px = 0.0f;
    float Py = 0.0f;
    float Pz = 0.0f;


    for (int i = 0; i < old_atom_accel_input.size(); i += 3)
    {
        float oldvx = old_velocity_input[i];
        float oldvy = old_velocity_input[i + 1];
        float oldvz = old_velocity_input[i + 2];

        float avgaccelx = 0.5 * (old_atom_accel_input[i] + new_atom_accel_input[i]) * dt;
        float avgaccely = 0.5 * (old_atom_accel_input[i + 1] + new_atom_accel_input[i + 1]) * dt;
        float avgaccelz = 0.5 * (old_atom_accel_input[i + 2] + new_atom_accel_input[i + 2]) * dt;

        float newvx = oldvx + avgaccelx;
        float newvy = oldvy + avgaccely;
        float newvz = oldvz + avgaccelz;

        // section to record kinetic energy and momentum
        KE_total += 0.5f * mass_input * ((newvx * newvx) + (newvy * newvy) + (newvz * newvz));
        Px += mass_input * newvx;
        Py += mass_input * newvy;
        Pz += mass_input * newvz;

        new_atom_velocity.push_back(newvx);
        new_atom_velocity.push_back(newvy);
        new_atom_velocity.push_back(newvz);
    }

    float momentum_mag = std::sqrt((Px * Px) + (Py * Py) + (Pz * Pz));
    momentum_magnitude.push_back(momentum_mag);
    kinetic_energy.push_back(KE_total);

    std::cout << "new_atom_velocity: " << std::endl;
    Print1DVec(new_atom_velocity);

    std::cout << "End of CalcVerletUpdateVelocity: " << std::endl;
    std::cout << std::endl;
    return new_atom_velocity;
}


/*
std::vector<float> Movement::ApplyPeriodicBoundaryConditions(std::vector<float> old_position_input,float box_length)
{
    std::cout << "Beginning of ApplyPeriodicBoundaryConditions: " << std::endl;

    new_pbc_atom_position.clear();

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
