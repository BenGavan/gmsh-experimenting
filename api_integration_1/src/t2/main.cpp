//
// Created by Ben Gavan on 14/07/2021.
//

#include <iostream>
#include "../../libs/gmsh.h"

using namespace std;

int main() {

    gmsh::initialize();

    gmsh::model::add("t2");

    // Copied from t1.cpp...
    double lc = 1e-2;
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(.1, 0, 0, lc, 2);
    gmsh::model::geo::addPoint(.1, .3, 0, lc, 3);
    gmsh::model::geo::addPoint(0, .3, 0, lc, 4);

    gmsh::model::geo::addLine(1, 2, 1);
    gmsh::model::geo::addLine(3, 2, 2);
    gmsh::model::geo::addLine(3, 4, 3);
    gmsh::model::geo::addLine(4, 1, 4);

    gmsh::model::geo::addCurveLoop({4, 1, -2, 3}, 1);

    gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::synchronize();

    gmsh::model::addPhysicalGroup(1, {1, 2, 4}, 5);
    int ps = gmsh::model::addPhysicalGroup(2, {1});
    gmsh::model::setPhysicalName(2, ps, "My surface");
    //


    return 0;
}

