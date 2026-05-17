#include "pdbreader.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // Required for fixed-width parsing (optional, can also use substr)
#include "../Utility/UtilityFunctions.h"

//https://www.google.com/search?q=c%2B%2B+code+to+read+a+protein+database+file&rlz=1C1VDKB_enUS1161US1161&oq=c%2B%2B+code+to+read+a+protein+database+file+&gs_lcrp=EgZjaHJvbWUyCQgAEEUYORifBTIHCAEQIRigATIHCAIQIRigATIHCAMQIRirAjIHCAQQIRirAjIHCAUQIRirAjIHCAYQIRifBTIHCAcQIRifBTIHCAgQIRifBTIHCAkQIRifBdIBCDY4MTVqMGo3qAIAsAIA&sourceid=chrome&ie=UTF-8

MyPdbInfo::MyPdbInfo()
{

}

//std::vector<float> pdbx; dont need this since im not giving it a number till in pdbinfo

// readPDB reads the pdb and creates a coordinate matrix
void MyPdbInfo::ReadPDB(const std::string& filename) {
    std::ifstream file(filename); //open the file "filename" for reading. file is the variable name
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    std::string line; // variable to store each line from the file 'file'
    while (std::getline(file, line)) { // reads each line from the file and stores it to print
        if (line.substr(0, 4) == "ATOM") { //HETATM
            // This is an ATOM record, parse relevant information

            std::string atom_index = line.substr(7,4);
            pdbatomindex.push_back(atom_index); // adds a new value to pdbatomindex while looping

            std::string atom_name = line.substr(12, 4); // Atom name, e.g., " CA "
            pdbatomname.push_back(atom_name); // adds a new value to pdbatomname while looping

            std::string res_name = line.substr(17, 3);  // Residue name, e.g., "ALA"
            pdbresname.push_back(res_name); // adds a new value to pdbresname while looping

            int res_seq = std::stoi(line.substr(22, 4)); // Residue sequence number
            // Coordinates
            float x = std::stof(line.substr(30, 8));
            float y = std::stof(line.substr(38, 8));
            float z = std::stof(line.substr(46, 8));
            x = x * std::pow(10,-1);        // need to convert from angstrom to nm
            y = y * std::pow(10,-1);        // need to convert from angstrom to nm
            z = z * std::pow(10,-1);        // need to convert from angstrom to nm

            pdbx.push_back(x);
            pdby.push_back(y);
            pdbz.push_back(z);

            // need to add a x1,y1,z1 column name to matrix eventually
            pdbmatrix.push_back({x,y,z});

            std::cout << "Atom: " << atom_name
                      << ", Index: " << atom_index
                      << ", Residue: " << res_name
                      << ", Seq: " << res_seq
                      << ", Coords: (" << x << ", " << y << ", " << z << ")" << std::endl;
        }
        // You can add more conditions for other record types (HEADER, REMARK, HETATM, etc.)
    }
    std::cout << std::endl;
    std::cout << "The name of the atoms are listed below in order by index: "<<std::endl;
    for (const auto& element:pdbatomname)
    {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    std::cout << "This is the pdbmatrix (x,y,z) in nm: "<< std::endl;
    for (const auto& row: pdbmatrix)
    {
        for (const auto& element: row)
        {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "This is pdbmatrix rows: "<< pdbmatrix.size() << std::endl;
    std::cout << "This is pdbmatrix columns: " << pdbmatrix[1].size() << std::endl;




    //////////////////////// printing 1d vector for pdbmatrix coordinates ///////////////////
    std::cout << std::endl;
    std::cout << "Vector elements for pdbx in format (atom1,atom2,atom3...): " << std::endl;
    // Use a range-based for loop to iterate through each element
    for (float element : pdbx) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    std::cout << "Vector elements for pdby in format (atom1,atom2,atom3...): " << std::endl;
    // Use a range-based for loop to iterate through each element
    for (float element : pdby) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    std::cout << "Vector elements for pdbz in format (atom1,atom2,atom3...): " << std::endl;
    // Use a range-based for loop to iterate through each element
    for (float element : pdbz) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
////////////////////////////////////////////////////////////////////////////////////////

    std::cout << std::endl;
    file.close();
    //readPDB("your_protein.pdb"); // Replace with your file path

}

//////////////////////////////////////////////////////////////
//going to create each section of the topology file then combine them into a final output
/////////////////////
////////////////////////////////////////////////////////////


/*

float MyPdbInfo::CalcDist3d(float x1, float y1, float z1, float x2, float y2, float z2)
{   // dist = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)

    float dx = x2-x1;
    float dy = y2-y1;
    float dz = z2-z1;

    return std::sqrt((dx*dx) + (dy*dy) + (dz*dz));


}
*/


void MyPdbInfo::MeasureBondDist(std::vector<std::vector<float>> pdbinput)
{
    // must be a symmetric matrix of form {{x1,y1,z1},{x2,y2,z2},...}
    std::cout << "Beginning MeasureBondDist Function: " << std::endl; // computes the distance hermitian matrix
    // have to do a double loop with pdbinput.size() so that every it loops i -> i+1, i-> i +2
    for (int i = 0; i < pdbinput.size();++i)
    {
        std::vector<float> result_matrix; // THIS IS THE RESULT_MATRIX
        //std::cout << "This is i: " << i << std::endl;

        for (int  j = 0; j < pdbinput.size(); ++j)
        {
            float result = CalcDist3d(
                pdbinput[i][0],pdbinput[i][1],pdbinput[i][2],
                pdbinput[j][0],pdbinput[j][1],pdbinput[j][2]);
            std::cout << "this is result: "<< result << std::endl;
            result_matrix.push_back(result);
            //std::cout << "this is j " << j << std::endl;

        }
        pdbdistmatrix.push_back(result_matrix);
        std::cout << std::endl;
    }

    // THis prints out the pdbdistmatrix
    std::cout <<"printing the pdbdistmaxtrix function section: "<< std::endl;

    for (const auto& element : pdbatomname)
    {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    //std::cout << std::endl;
    for (const auto& row : pdbdistmatrix) // prints the 2d pdbdistmatrix
    {
        for (float element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

}

std::vector<std::vector<float>> MyPdbInfo::Pdb2Gmx(std::vector<std::vector<float>> pdbinput)
{
    std::cout << "Beginning of Pdb2Gmx function: " << std::endl;
    std::cout << "Converts the pdbinput from {{x1,y1,z1,},{x2,y2,z2}} to {{x1,y1,z1,x2,y2,z2,x3....} " << std::endl;

    std::cout << "The pdbinput.size() is " << pdbinput.size() << std::endl;      // rows
    std::cout << "The pdbinput[0].size() is "<< pdbinput[0].size() << std::endl; // columns


    pdb2gmxmatrix={pdbinput[0]};
    for (int i = 1; i < pdbinput.size(); ++i) // starts at i=1 because .inser() has to have some starting value
    {
        pdb2gmxmatrix[0].insert(pdb2gmxmatrix[0].end(),pdbinput[i].begin(),pdbinput[i].end());
    }

    std::cout << "The pdb2gmx matrix is " << std::endl;
    std::cout << "The pdb2gmx matrix.size() is: " << pdb2gmxmatrix.size() << std::endl;
    std::cout << "The pdb2gmx matrix[0].size() is: " << pdb2gmxmatrix[0].size() << std::endl;
    for (const auto& row : pdb2gmxmatrix) // prints the 2d pdb2gmxmatrix
    {
        for (float element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    return pdb2gmxmatrix;
}

std::vector<std::vector<float>> MyPdbInfo::GmxMeasureDist(std::vector<std::vector<float>> pdb2gmxmatrix_input)
{
    std::cout << "Beginning of GmxMeasureDist Function: " << std::endl;
    std::cout << std::endl;

    const std::vector<float>& last_row = pdb2gmxmatrix_input.back();
    for (int i=0;i < last_row.size();i+=3)
    {
        std::vector<float> result_matrix;
        for (int j=0; j< last_row.size();j+=3)
        {
            float result = CalcDist3d(
                last_row[i],last_row[i+1],last_row[i+2],
                last_row[j],last_row[j+1],last_row[j+2]);
            std::cout << "This is result from GmxMeasureBondDist: " << result << std::endl;
            result_matrix.push_back(result);
        }
        pdb2gmx_distmatrix.push_back(result_matrix);
        std::cout << std::endl;
    }

    std::cout << "Printing the pdb2gmx_distmatrix: " << std::endl;
    for (const auto& element: pdbatomname) // prints the atom names first
    {
        std::cout << element << " ";
    }
    std::cout << std::endl;

    for (const auto& row : pdb2gmx_distmatrix)
    {
        for (float element : row)
        {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    } //prints the pdb2gmx_distmatrix
    std::cout << "End of pdbreader.cpp:: GmxMeasureDist" << std::endl;
    std::cout << std::endl;


    return pdb2gmx_distmatrix;

}