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

    // create a new view
    int t1 = gmsh::view::add("A list-based view");

    // create a triangle
    vector<double> triangle1 = {
            0., 1., 1., // x coordinates of the 3 triangle nodes
            0., 0., 1., // y coordinates of the 3 triangle nodes
            0., 0., 0., // z coordinates of the 3 triangle nodes
    };

    // add values at each node
    triangle1.insert(triangle1.end(), {10, 11, 12});

    // concatenate data
    vector<double> triangles(triangle1);
//    triangles.insert(triangles.end(), );
    gmsh::view::addListData(t1, "ST", 1, triangles);

    gmsh::view::write(t1, "x3.pos");

    gmsh::finalize();
}
