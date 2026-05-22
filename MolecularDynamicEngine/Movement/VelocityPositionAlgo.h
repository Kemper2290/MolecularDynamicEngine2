#ifndef VELOCITYPOSITIONALGO_H
#define VELOCITYPOSITIONALGO_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>


class Movement
{
public:

    // Member Function
    static std::vector<float> CalcAcceleration(float mass_input,
                                        std::vector<float> atom_force_input);



    static std::vector<float> CalcVerletPosition(std::vector<std::vector<float>> coord_input
                                                            ,std::vector<float> atom_accel_input,
                                                            std::vector<float> velocity_input,
                                                            float dt,
                                                            bool pbc,
                                                            float box_length_input);

    static std::vector<float> PBCCalcVerletPosition(std::vector<std::vector<float>> position_inp,
                                                        std::vector<float> velocity_inp,
                                                        std::vector<float> acceleration_inp,
                                                        float dt,
                                                        float box_length,
                                                        bool pbc,
                                                        bool rbc);


    static std::vector<float> CalcVerletUpdateVelocity(std::vector<float> old_velocity_input,
                                                        std::vector<float> old_atom_accel_input,
                                                        std::vector<float> new_atom_accel_input,
                                                        float dt,
                                                        float mass_input);


    static void ApplyReflectiveWalls(std::vector<float>& position, std::vector<float>& velocity,float boxsize);

    static void CenterCoordinatesInBox(std::vector<float>& coords, float box_length);

    // variables


    //static std::vector<float> force;
    static std::vector<float> atom_accel;

    static std::vector<float> new_atom_position;
    static std::vector<float> new_atom_velocity;
    //static std::vector<float> new_pbc_atom_position;
    static std::vector<float> kinetic_energy;
    static std::vector<float> momentum_magnitude;



private:




};










#endif