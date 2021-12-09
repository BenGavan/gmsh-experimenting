//
// Created by Ben Gavan on 18/07/2021.
//

#include <iostream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;

int main()
{
    int t1 = gmsh::view::add("A list-based view");

    cout << "t1: " << t1 << '\n';

    std::vector<double> triangle1 = {
            0., 1., 1., // x coordinates of the 3 triangle nodes
            0., 0., 1., // y coordinates of the 3 triangle nodes
            0., 0., 0.}; // z coordinates of the 3 triangle nodes
    std::vector<double> triangle2 = {0., 1., 0., 0., 1., 1., 0., 0., 0.};


    for(int step = 0; step < 10; step++) {
        triangle1.insert(triangle1.end(), {10., 11. - step, 12.});
        triangle2.insert(triangle2.end(), {11., 12., 13. + step});
    }

    // List-based data is just added by concatenating the data for all the
    // triangles:
    std::vector<double> triangles(triangle1);
    triangles.insert(triangles.end(), triangle2.begin(), triangle2.end());
    gmsh::view::addListData(t1, "ST", 2, triangles);

    gmsh::model::getBoundary()
    return 0;
}

