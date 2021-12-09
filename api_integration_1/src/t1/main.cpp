//
// Created by Ben Gavan on 14/07/2021.
//

/*
 * Makes a square then meshes it with default parameters
 */

#include <iostream>
#include "../../libs/gmsh.h"

int main() {

    gmsh::initialize();

    gmsh::model::add("t1");

    double lc = 1e-2;
    // points
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(.1, 0, 0, lc, 2);
    gmsh::model::geo::addPoint(.1, .3, 0, lc, 3);
    gmsh::model::geo::addPoint(0, .3, 0, lc, 4);

    // lines (from points)
    gmsh::model::geo::addLine(1, 2, 1);
    gmsh::model::geo::addLine(3, 2, 2);
    gmsh::model::geo::addLine(3, 4, 3);
    gmsh::model::geo::addLine(4, 1, 4);

    // loops (from lines) (to make surfaces out of)
    gmsh::model::geo::addCurveLoop({4, 1, -2, 3}, 1);

    // surfaces (from loops)
    gmsh::model::geo::addPlaneSurface({1}, 1);

    // needs to be called before used outside of CAD kernel (e.g. meshing)
    gmsh::model::geo::synchronize();

    // create physical group
//    gmsh::model::addPhysicalGroup(1, {1, 2, 4}, 5);
//    int ps = gmsh::model::addPhysicalGroup(2, {1});
//    gmsh::model::setPhysicalName(2, ps, "My Surface");

    // generate mesh
    gmsh::model::mesh::generate(2);

    gmsh::write("../../out/t1.msh");

    return 0;
}

