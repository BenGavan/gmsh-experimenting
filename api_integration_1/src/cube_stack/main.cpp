//
// Created by Ben Gavan on 20/07/2021.
//

#include <iostream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;


void make_cube(const double x_shift)
{
//    string name = "t1_cube_";
//    name.append(to_string(x_shift));
//    gmsh::model::add(name);

    double lc = 1e-2;

    // Make cube of points
    int p1 = gmsh::model::geo::addPoint(0+x_shift, 0, 0, lc);
    int p2 = gmsh::model::geo::addPoint(0+x_shift, 1, 0, lc);
    int p3 = gmsh::model::geo::addPoint(0+x_shift, 0, 1, lc);
    int p4 = gmsh::model::geo::addPoint(0+x_shift, 1, 1, lc);

    int p5 = gmsh::model::geo::addPoint(1+x_shift, 0, 0, lc);
    int p6 = gmsh::model::geo::addPoint(1+x_shift, 1, 0, lc);
    int p7 = gmsh::model::geo::addPoint(1+x_shift, 0, 1, lc);
    int p8 = gmsh::model::geo::addPoint(1+x_shift, 1, 1, lc);

    // Make lines to form cube
    int l1 = gmsh::model::geo::addLine(p1, p2);
    int l2 = gmsh::model::geo::addLine(p1,p3);
    int l3 = gmsh::model::geo::addLine(p1,p5);
    int l4 = gmsh::model::geo::addLine(p2,p6);
    int l5 = gmsh::model::geo::addLine(p2,p4);
    int l6 = gmsh::model::geo::addLine(p6,p5);
    int l7 = gmsh::model::geo::addLine(p6,p8);
    int l8 = gmsh::model::geo::addLine(p8,p4);
    int l9 = gmsh::model::geo::addLine(p8,p7);
    int l10 = gmsh::model::geo::addLine(p4, p3);
    int l11 = gmsh::model::geo::addLine(p7, p3);
    int l12 = gmsh::model::geo::addLine(p5, p7);

    // Make curve loops to be able to make surfaces from them
    int c1 = gmsh::model::geo::addCurveLoop({l1, l5, l10, -l2});
    int c2 = gmsh::model::geo::addCurveLoop({l4, l7, l8, -l5});
    int c3 = gmsh::model::geo::addCurveLoop({l1, l4, l6, -l3});
    int c4 = gmsh::model::geo::addCurveLoop({l3, l12, l11, -l2});
    int c5 = gmsh::model::geo::addCurveLoop({l6, l12, -l9, -l7});
    int c6 = gmsh::model::geo::addCurveLoop({l8, l10, -l11, -l9});

    // Make surfaces from curve loops
    int ps1 = gmsh::model::geo::addPlaneSurface({c1});
    int ps2 = gmsh::model::geo::addPlaneSurface({c2});
    int ps3 = gmsh::model::geo::addPlaneSurface({c3});
    int ps4 = gmsh::model::geo::addPlaneSurface({c4});
    int ps5 = gmsh::model::geo::addPlaneSurface({c5});
    int ps6 = gmsh::model::geo::addPlaneSurface({c6});

    // Make surface loop from plane surfaces
    int sl1 = gmsh::model::geo::addSurfaceLoop({ps1,ps2,ps3,ps4,ps5,ps6});

    // Make volume from surface loop
    gmsh::model::geo::addVolume({sl1});
}

int main()
{
    gmsh::initialize();

    try {
        gmsh::merge("bg.pos");
    } catch(...) {
        gmsh::logger::write("Could not load background mesh: bye!", "error");
        gmsh::finalize();
        return 0;
    }

    gmsh::model::add("cube");
    for (double i=0; i<40; i++) {
        make_cube(i);
    }

    gmsh::model::geo::synchronize();

    // Add the post-processing view as a new size field:
    int bg_field = gmsh::model::mesh::field::add("PostView");
    gmsh::model::mesh::field::setNumber(bg_field, "ViewIndex", 0);

    gmsh::model::mesh::field::setAsBackgroundMesh(bg_field);

    // In order to compute the mesh sizes from the background mesh only, and
    // disregard any other size constraints, one can set:
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0); // 0 => mesh generation not based on Mesh.x
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    gmsh::model::mesh::generate(3);

    gmsh::write("../../out/cube_stack.msh");

    // Always finish with finalize
    gmsh::finalize();


    return 0;
}

