#ifndef UTILITYFUNCTIONS_h
#define UTILITYFUNCTIONS_h

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <variant>
#include <chrono>
#include <random>
#include <cmath>


template <typename T>
void Print1DVec(T vec);

template <typename T>
void Print2DVec(T vec);

float CalcDist3d(float x1, float y1, float z1, float x2, float y2, float z2);

float PBCCalcDist3d(float x1, float y1, float z1, float x2, float y2, float z2, float box_length);

std::vector<float> Gen_Velocity(std::vector<std::vector<float>> pdb2gmx_coord_input,
                                                float mass_input,
                                                float temp_k);

void Plot_Trajectory(std::vector<std::vector<float>> trajectory_input,float delay,float xlim,float ylim);

void Plot2D(std::vector<std::vector<float>> trajectory_input,float delay, float xlim,float ylim);

//std::vector<float> ApplyPeriodicBoundaryConditions(std::vector<float> old_position_input,float box_lenght);
//std::vector<float> new_pbc_atom_position;

std::vector<std::vector<float>> CreateAtomCluster(float space, float gridsize = 3);

std::vector<float> MIC_2DTrajectory(std::vector<float>& position_inp, float box_size_inp);

#endif