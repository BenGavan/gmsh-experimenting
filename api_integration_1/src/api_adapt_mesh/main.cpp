//
// Created by Ben Gavan on 20/07/2021.
//

#include <iostream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;

int main()
{
    double lc = 0.02;
    int N = 10000;

    gmsh::initialize();

    // create a geometrical gmsh.model
    gmsh::model::add("square");
    int square = gmsh::model::occ::addRectangle(0, 0, 0, 1, 1);
    gmsh::model::occ::synchronize();

    // create intial uniform mesh
    gmsh::vectorpair in = {pair<int,int>(2, square)};
    gmsh::vectorpair pnts;

    gmsh::model::getBoundary(in, pnts, true, true, true);
    gmsh::model::mesh::setSize(pnts, lc);
    gmsh::model::mesh::generate(2);
    gmsh::write("../../out/api_adapt_mesh.msh");

    gmsh::finalize();

    return 0;
}
