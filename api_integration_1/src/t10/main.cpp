//
// Created by Ben Gavan on 16/07/2021.
//

#include <set>
#include <sstream>
#include <iostream>
//#include <gmsh.h>
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

    gmsh::model::mesh::field::add("Distance", 1);
    gmsh::model::mesh::field::setNumbers(1, "PointsList", {5});
    gmsh::model::mesh::field::setNumbers(1, "CurvesList", {2});
    gmsh::model::mesh::field::setNumber(1, "NumPointsPerCurve", 100);

    // We then define a `Threshold' field, which uses the return value of the
    // `Distance' field 1 in order to define a simple change in element size
    // depending on the computed distances
    //
    // SizeMax -                     /------------------
    //                              /
    //                             /
    //                            /
    // SizeMin -o----------------/
    //          |                |    |
    //        Point         DistMin  DistMax
    gmsh::model::mesh::field::add("Threshold", 2);
    gmsh::model::mesh::field::setNumber(2, "InField", 1);
    gmsh::model::mesh::field::setNumber(2, "SizeMin", lc / 30);
    gmsh::model::mesh::field::setNumber(2, "SizeMax", lc);
    gmsh::model::mesh::field::setNumber(2, "DistMin", 0.15);
    gmsh::model::mesh::field::setNumber(2, "DistMax", 0.5);

//    // Say we want to modulate the mesh element sizes using a mathematical
//    // function of the spatial coordinates. We can do this with the MathEval
//    // field:
    gmsh::model::mesh::field::add("MathEval", 3);
    gmsh::model::mesh::field::setString(3, "F", "Cos(4*3.14*x) * Sin(4*3.14*y) / 10 + 0.101");

    // We could also combine MathEval with values coming from other fields. For
    // example, let's define a `Distance' field around point 1
    gmsh::model::mesh::field::add("Distance", 4);
    gmsh::model::mesh::field::setNumbers(4, "PointsList", {1, 5, 4});

    // We can then create a `MathEval' field with a function that depends on the
    // return value of the `Distance' field 4, i.e., depending on the distance to
    // point 1 (here using a cubic law, with minimum element size = lc / 100)
    gmsh::model::mesh::field::add("MathEval", 5);
    std::stringstream stream;
    stream << "F4^3 + " << lc / 100;
    gmsh::model::mesh::field::setString(5, "F", stream.str());

//    // We could also use a `Box' field to impose a step change in element sizes
//    // inside a box
    gmsh::model::mesh::field::add("Box", 6);
    gmsh::model::mesh::field::setNumber(6, "VIn", lc / 15);
    gmsh::model::mesh::field::setNumber(6, "VOut", lc);
    gmsh::model::mesh::field::setNumber(6, "XMin", 0.3);
    gmsh::model::mesh::field::setNumber(6, "XMax", 0.6);
    gmsh::model::mesh::field::setNumber(6, "YMin", 0.3);
    gmsh::model::mesh::field::setNumber(6, "YMax", 0.6);
//
////    // Many other types of fields are available: see the reference manual for a
////    // complete list. You can also create fields directly in the graphical user
////    // interface by selecting `Define->Size fields' in the `Mesh' module.
////
////    // Finally, let's use the minimum of all the fields as the background mesh
//    // field:
    gmsh::model::mesh::field::add("Min", 7);
    gmsh::model::mesh::field::setNumbers(7, "FieldsList", {2, 3, 5, 6});
//
    gmsh::model::mesh::field::setAsBackgroundMesh(7);

//    // The API also allows to set a global mesh size callback, which is called
//    // each time the mesh size is queried
//    auto meshSizeCallback = [](int dim, int tag, double x, double y, double z) {
//        return 0.02 * x * y + 0.01;
//        return 0.02 * x + 0.01;
//    };
//    gmsh::model::mesh::setSizeCallback(meshSizeCallback);

    // to just using the background mesh
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);


    gmsh::model::mesh::generate(2);
    gmsh::write("../../out/t10.msh");

    gmsh::finalize();
}

