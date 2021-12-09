//
// Created by Ben Gavan on 24/07/2021.
//

#include <iostream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;

int main()
{
    gmsh::initialize();

    gmsh::model::add("t10");

    // create a simple rectangular geometry:
    double lc = 2;
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

    // create a new view
    int view_tag = gmsh::view::add("A list-based view");

    cout << view_tag << endl;

    // create a triangle
    vector<double> triangle1 = {
            0., 1., 1., // x coordinates of the 3 triangle nodes
            0., 0., 1., // y coordinates of the 3 triangle nodes
            0., 0., 0., // z coordinates of the 3 triangle nodes
    };

    vector<double> triangle2 = {
            0., 1., 1., // x coordinates of the 3 triangle nodes
            0., 0., 1., // y coordinates of the 3 triangle nodes
            0., 0., 0., // z coordinates of the 3 triangle nodes
    };

    for (double x : triangle1)
    {
        cout << x << ", ";
    }
    cout << endl;

    // add values at each node
    triangle1.insert(triangle1.end(), {.01, .0005, 1.2});

    for (double x : triangle1)
    {
        cout << x << ", ";
    }
    cout << endl;

    // concatenate data
    vector<double> triangles(triangle1);
//    triangles.insert(triangles.end(), {1});
//    triangles.insert(triangles.end(), triangle2.begin(),  triangle2.end());
    gmsh::view::addListData(view_tag, "ST", 1, triangles);

    gmsh::view::write(view_tag, "out.pos");


    int post_tag = gmsh::model::mesh::field::add("PostView");
    gmsh::model::mesh::field::setNumber(post_tag, "ViewTag", view_tag);

    gmsh::model::mesh::field::add("Box", 2);
    gmsh::model::mesh::field::setNumber(2, "VIn", lc / 150);
    gmsh::model::mesh::field::setNumber(2, "VOut", lc);
    gmsh::model::mesh::field::setNumber(2, "XMin", 0.);
    gmsh::model::mesh::field::setNumber(2, "XMax", 0.25);
    gmsh::model::mesh::field::setNumber(2, "YMin", 0.75);
    gmsh::model::mesh::field::setNumber(2, "YMax", 1.);

    gmsh::model::mesh::field::add("Min", 7);
    gmsh::model::mesh::field::setNumbers(7, "FieldsList", {(double)post_tag,2});

    gmsh::model::mesh::field::setAsBackgroundMesh(7);
//
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

//    gmsh::model::mesh::field::
    gmsh::model::mesh::generate(2);
    gmsh::write("out.msh");

    gmsh::finalize();
}
