//
// Created by Ben Gavan on 13/07/2021.
//

/*
 * Creates cubes from individual points
 * Then meshes the surface of cube with default parameters.
 */

#include <iostream>
#include "../../libs/gmsh.h"

using namespace std;

void make_cube()
{
    gmsh::model::add("t1_cube");

    double lc = 1e-2;

    // Make cube of points
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(0, 1, 0, lc, 2);
    gmsh::model::geo::addPoint(0, 0, 1, lc, 3);
    gmsh::model::geo::addPoint(0, 1, 1, lc, 4);

    gmsh::model::geo::addPoint(1, 0, 0, lc, 5);
    gmsh::model::geo::addPoint(1, 1, 0, lc, 6);
    gmsh::model::geo::addPoint(1, 0, 1, lc, 7);
    gmsh::model::geo::addPoint(1, 1, 1, lc, 8);

    // Make lines to form cube
    gmsh::model::geo::addLine(1, 2, 1);
    gmsh::model::geo::addLine(1,3,2);
    gmsh::model::geo::addLine(1,5,3);
    gmsh::model::geo::addLine(2,6,4);
    gmsh::model::geo::addLine(2,4,5);
    gmsh::model::geo::addLine(6,5,6);
    gmsh::model::geo::addLine(6,8,7);
    gmsh::model::geo::addLine(8,4,8);
    gmsh::model::geo::addLine(8,7,9);
    gmsh::model::geo::addLine(4, 3,  10);
    gmsh::model::geo::addLine(7, 3, 11);
    gmsh::model::geo::addLine(5, 7, 12);

    // Make curve loops to be able to make surfaces from them
    gmsh::model::geo::addCurveLoop({1, 5, 10, -2}, 1);
    gmsh::model::geo::addCurveLoop({4, 7, 8, -5}, 2);
    gmsh::model::geo::addCurveLoop({1, 4, 6, -3}, 3);
    gmsh::model::geo::addCurveLoop({3, 12, 11, -2}, 4);
    gmsh::model::geo::addCurveLoop({6, 12, -9, -7}, 5);
    gmsh::model::geo::addCurveLoop({8, 10, -11, -9}, 6);

    // Make surfaces from curve loops
    gmsh::model::geo::addPlaneSurface({1}, 1);
    gmsh::model::geo::addPlaneSurface({2}, 2);
    gmsh::model::geo::addPlaneSurface({3}, 3);
    gmsh::model::geo::addPlaneSurface({4}, 4);
    gmsh::model::geo::addPlaneSurface({5}, 5);
    gmsh::model::geo::addPlaneSurface({6}, 6);

    // Make surface loop from plane surfaces
    gmsh::model::geo::addSurfaceLoop({1,2,3,4,5,6}, 1);

    // Make volume from surface loop
    gmsh::model::geo::addVolume({1}, 1);

    gmsh::model::geo::synchronize();
}

int main() {
    gmsh::initialize();

    make_cube();

    auto meshSizeCallback = [](int dim, int tag, double x, double y, double z) {
//        return pow(0.5, 3) /  8.485;
        return .25;
    };
    gmsh::model::mesh::setSizeCallback(meshSizeCallback);

    // to just using the background mesh
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);


    gmsh::model::mesh::generate(3);

    gmsh::write("../../out/cube_withcallback_0.25.msh");

    // Always finish with finalize
    gmsh::finalize();

    return 0;
}