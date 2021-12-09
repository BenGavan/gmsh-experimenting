//
// Created by Ben Gavan on 18/07/2021.
//

#include <iostream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;

int main()
{
    gmsh::initialize();

    // open msh
    string filename = "t7.msh"; // .msh
    gmsh::open(filename);

    std::string name;
    gmsh::model::getCurrent(name);
    std::cout << "Model " << name << " (" << gmsh::model::getDimension() << "D)\n";

    // close msh
    gmsh::clear();

    // open background msh (.pos)
    filename = "t7_bgmesh.pos"; // background mesh
    gmsh::open(filename);
    gmsh::model::getCurrent(name);
    cout << "Model 2: " << name << " (" << gmsh::model::getDimension() << "D)\n";

    std::vector<std::pair<int, int> > entities;
    gmsh::model::getEntities(entities);

    cout << entities.size() << endl;

    gmsh::clear();

    gmsh::finalize();
    return 0;
}


