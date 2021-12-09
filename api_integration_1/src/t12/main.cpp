//
// Created by Ben Gavan on 16/07/2021.
//

#include <set>
#include <sstream>
#include <iostream>
#include "../../libs/gmsh.h"

int main()
{
    gmsh::initialize();

    gmsh::model::add("t12");

    double lc = 0.1;

    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(1, 0, 0, lc, 2);
    gmsh::model::geo::addPoint(1, 1, 0.5, lc, 3);
    gmsh::model::geo::addPoint(0, 1, 0.4, lc, 4);
    gmsh::model::geo::addPoint(0.3, 0.2, 0, lc, 5);
    gmsh::model::geo::addPoint(0, 0.01, 0.01, lc, 6);
    gmsh::model::geo::addPoint(0, 0.02, 0.02, lc, 7);
    gmsh::model::geo::addPoint(1, 0.05, 0.02, lc, 8);
    gmsh::model::geo::addPoint(1, 0.32, 0.02, lc, 9);

    gmsh::model::geo::addLine(1, 2, 1);
    gmsh::model::geo::addLine(2, 8, 2);
    gmsh::model::geo::addLine(8, 9, 3);
    gmsh::model::geo::addLine(9, 3, 4);
    gmsh::model::geo::addLine(3, 4, 5);
    gmsh::model::geo::addLine(4, 7, 6);
    gmsh::model::geo::addLine(7, 6, 7);
    gmsh::model::geo::addLine(6, 1, 8);
    gmsh::model::geo::addSpline({7, 5, 9}, 9);
    gmsh::model::geo::addLine(6, 8, 10);

    gmsh::model::geo::addCurveLoop({5, 6, 9, 4}, 11);
    gmsh::model::geo::addSurfaceFilling({11}, 1);

    gmsh::model::geo::addCurveLoop({-9, 3, 10, 7}, 13);
    gmsh::model::geo::addSurfaceFilling({13}, 5);

    gmsh::model::geo::addCurveLoop({-10, 2, 1, 8}, 15);
    gmsh::model::geo::addSurfaceFilling({15}, 10);

    gmsh::model::geo::synchronize();

    // Treat curves 2, 3 and 4 as a single curve when meshing (i.e. mesh across
    // points 6 and 7)
    gmsh::model::mesh::setCompound(1, {2, 3, 4});

    // Idem with curves 6, 7 and 8
    gmsh::model::mesh::setCompound(1, {6, 7, 8});

    // Treat surfaces 1, 5 and 10 as a single surface when meshing (i.e. mesh
    // across curves 9 and 10)
    gmsh::model::mesh::setCompound(2, {1, 5, 10});

    gmsh::model::mesh::generate(2);
    gmsh::write("../../out/t12.msh");

    gmsh::finalize();

    return 0;
}