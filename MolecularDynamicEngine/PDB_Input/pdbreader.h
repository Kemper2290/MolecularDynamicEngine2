#ifndef PDBBREADER_H
#define PDBBREADER_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

// function declaration

class MyPdbInfo
{
public:
    // public members ( accessible from outside the class
    MyPdbInfo(); // constructor declaration

    // Methods
    void ReadPDB(const std::string& filename); //public function declaration
    void MeasureBondDist(std::vector<std::vector<float>> pdbinput);
    //float CalcDist3d(float x1,float y1, float z1, float x2, float y2, float z2);
    std::vector<std::vector<float>> Pdb2Gmx(std::vector<std::vector<float>> pdbinput);
    std::vector<std::vector<float>> GmxMeasureDist(std::vector<std::vector<float>> pdb2gmxmatrix_input);

    // Variables
    std::vector<std::vector<float>> pdbmatrix; // prints the x,y,z matrix for atom coordinates {{x1,y1,z1},{x2,y2,z2}...}
    std::vector<std::string> pdbatomindex; // stores the atom index from pdb input
    std::vector<std::string> pdbatomname; // stores atom names from pdb input
    std::vector<std::string> pdbresname; // stores residue names from pdb input
    std::vector<std::vector<float>> pdbbonds;
    std::vector<std::vector<float>> pdbdistmatrix;
    std::vector<float> pdbx;
    std::vector<float> pdby;
    std::vector<float> pdbz;
    std::vector<std::vector<float>> pdb2gmxmatrix;
    std::vector<std::vector<float>> pdb2gmx_distmatrix;

private:
    // private members( accesible only within the class


};




#endif
