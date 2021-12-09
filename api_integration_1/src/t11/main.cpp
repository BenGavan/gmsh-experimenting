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

    gmsh::model::add("t11");

    // We have seen in tutorials `t3.cpp' and `t6.cpp' that extruded and
    // transfinite meshes can be "recombined" into quads, prisms or
    // hexahedra. Unstructured meshes can be recombined in the same way. Let's
    // define a simple geometry with an analytical mesh size field:

    int p1 = gmsh::model::geo::addPoint(-1.25, -.5, 0);
    int p2 = gmsh::model::geo::addPoint(1.25, -.5, 0);
    int p3 = gmsh::model::geo::addPoint(1.25, 1.25, 0);
    int p4 = gmsh::model::geo::addPoint(-1.25, 1.25, 0);

    int l1 = gmsh::model::geo::addLine(p1, p2);
    int l2 = gmsh::model::geo::addLine(p2, p3);
    int l3 = gmsh::model::geo::addLine(p3, p4);
    int l4 = gmsh::model::geo::addLine(p4, p1);

    int cl = gmsh::model::geo::addCurveLoop({l1, l2, l3, l4});
    int pl = gmsh::model::geo::addPlaneSurface({cl});

    gmsh::model::geo::synchronize();

    gmsh::model::mesh::field::add("MathEval", 1);
    gmsh::model::mesh::field::setString(
            1, "F", "0.01*(1.0+30.*(y-x*x)*(y-x*x) + (1-x)*(1-x))");
    gmsh::model::mesh::field::setAsBackgroundMesh(1);

    // To generate quadrangles instead of triangles, we can simply add
//    gmsh::model::mesh::setRecombine(2, pl);

//    gmsh::option::setNumber("Mesh.Algorithm", 8); //  for a preference of right-angles.

    gmsh::model::mesh::generate(2);

    gmsh::model::mesh::recombine();
    gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);
    gmsh::model::mesh::refine();

    gmsh::write("../../out/t11.msh");

    gmsh::finalize();
}

