//
// Created by Ben Gavan on 22/07/2021.
//

#include <iostream>
#include <fstream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;

vector<size_t> getTetrahedronMeshElementTags()
{
    vector<size_t> elementTags;
    vector<size_t> nodeTags;
    gmsh::model::mesh::getElementsByType(4, elementTags, nodeTags); // 4 = 4-node tetrahedron
    return elementTags;
}

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

int main()
{
    gmsh::initialize();

    for (unsigned i=0; i<10; i++) {
        // setup unit cube
        make_cube();

        // mesh
        double size = 1 / (pow(2, i));

        auto meshSizeCallback = [=](int dim, int tag, double x, double y, double z) {
            return size;
        };
        gmsh::model::mesh::setSizeCallback(meshSizeCallback);

        // to just using the background mesh
        gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
        gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
        gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

        gmsh::model::mesh::generate(3);

        // count tetrahedron elements
        vector<size_t> tags = getTetrahedronMeshElementTags();
        int num_elements = tags.size();

        // write to file
        cout << "size: " << 1 / (pow(2, i)) << ", #elements: " << num_elements << endl;
        ofstream f;
        f.open("out.txt", ios_base::app);
        f << size << "," << num_elements << '\n';
        f.close();
    }

    gmsh::finalize();
    return 0;
}
