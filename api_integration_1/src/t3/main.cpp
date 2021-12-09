//
// Created by Ben Gavan on 14/07/2021.
//

#include <iostream>
#include "../../libs/gmsh.h"

using namespace std;

int main() {
    gmsh::initialize();

    auto createGeometryAndMesh = []() {
        // Clear all models and create a new one
        gmsh::clear();
        gmsh::model::add("t3");

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


        double h = 0.1;
        std::vector <std::pair<int, int>> ov;
        gmsh::model::geo::extrude({{2, 1}}, 0, 0, h, ov, {8, 2}, {0.5, 1});

        // The extrusion can also be performed with a rotation instead of a
        // translation, and the resulting mesh can be recombined into prisms (we use
        // only one layer here, with 7 subdivisions). All rotations are specified by
        // an an axis point (-0.1, 0, 0.1), an axis direction (0, 1, 0), and a
        // rotation angle (-Pi/2):
        gmsh::model::geo::revolve({{2, 28}}, -0.1, 0, 0.1, 0, 1, 0, -M_PI / 2, ov,
                                  {7});

        std::vector<double> angle;
        gmsh::onelab::getNumber("Parameters/Twisting angle", angle);
        gmsh::model::geo::twist({{2, 50}}, 0, 0.15, 0.25, -2 * h, 0, 0, 1, 0, 0,
                                angle[0] * M_PI / 180., ov, {10}, {}, true);

        gmsh::model::geo::synchronize();

        // All the extrusion functions return a vector of extruded entities: the
        // "top" of the extruded surface (in `ov[0]'), the newly created volume (in
        // `ov[1]') and the tags of the lateral surfaces (in `ov[2]', `ov[3]', ...).

        // We can then define a new physical volume (with tag 101) to group all the
        // elementary volumes:
        gmsh::model::addPhysicalGroup(3, {1, 2, ov[1].second}, 101);

        gmsh::model::mesh::generate(3);
//        gmsh::write("../../out/t3.msh");
        gmsh::write("t3.msh");
    };


    gmsh::option::setNumber("Geometry.PointNumbers", 1);
    gmsh::option::setColor("Geometry.Points", 255, 165, 0);
    gmsh::option::setColor("General.Text", 255, 255, 255);
    gmsh::option::setColor("Mesh.Points", 255, 0, 0);

    int r, g, b, a;
    gmsh::option::getColor("Geometry.Points", r, g, b, a);
    gmsh::option::setColor("Geometry.Surfaces", r, g, b, a);

    // Twist parameters
    gmsh::onelab::set(R"( [
  {
    "type":"number",
    "name":"Parameters/Twisting angle",
    "values":[90],
    "min":0,
    "max":120,
    "step":1
  }
  ] )");

    createGeometryAndMesh();
}

