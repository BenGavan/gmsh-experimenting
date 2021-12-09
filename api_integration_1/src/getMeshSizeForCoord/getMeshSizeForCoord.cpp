//
// Created by Ben Gavan on 16/07/2021.
//



double getSizeForCoord(x double, y double, z, double) {
    return x * y * z;
}

//
// Created by Ben Gavan on 16/07/2021.
//

#include <set>
#include <sstream>
#include <iostream>
#include "../../libs/gmsh.h"

using namespace std;

int main()
{
    gmsh::initialize();

    gmsh::model::add("t10");

    // create a simple rectangular geometry:
    double lc = .15;
    gmsh::model::geo::addPoint(0.0, 0.0, 0, lc, 1);
    gmsh::model::geo::addPoint(1, 0.0, 0, lc, 2);
    gmsh::model::geo::addPoint(1, 1, 0, lc, 3);
    gmsh::model::geo::addPoint(0, 1, 0, lc, 4);
    gmsh::model::geo::addPoint(0.2, .5, 0, lc, 5);

    gmsh::model::geo::addLine(1, 2, 1);
    gmsh::model::geo::addLine(2, 3, 2);
    gmsh::model::geo::addLine(3, 4, 3);
    gmsh::model::geo::addLine(4, 1, 4);

    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 5);
    gmsh::model::geo::addPlaneSurface({5}, 6);

    gmsh::model::geo::synchronize();


    // The API also allows to set a global mesh size callback, which is called
    // each time the mesh size is queried
    auto meshSizeCallback = [](int dim, int tag, double x, double y, double z) {
        return getSizeForCoord(x, y, x);
    };
    gmsh::model::mesh::setSizeCallback(meshSizeCallback);

    // so just using the background mesh
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);


    gmsh::model::mesh::generate(2);
    gmsh::write("../../out/t10.msh");

    gmsh::finalize();
}


