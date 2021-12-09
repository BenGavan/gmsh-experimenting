//
// Created by Ben Gavan on 23/07/2021.
//

#include <iostream>
#include <fstream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;


// Node //
class Node {
public:
    size_t tag;
    std::vector<double> coord;
    std::vector<double> parametricCoord;
    double target_size;

    Node(std::size_t, std::vector<double>, std::vector<double>);
    void printCoords();
};

Node::Node(std::size_t t, std::vector<double> c, std::vector<double> pc)
{
    tag = t;
    coord = c;
    parametricCoord = pc;
}

void Node::printCoords()
{
    cout << "(" << coord[0] << ", " << coord[1] << ", " << coord[2] << ")";
}
// //

// Element //
class Element
{
public:
    size_t tag;
    size_t type;
    vector<Node> nodes;

    Element(size_t tag, size_t type, vector<Node> nodes);
};

Element::Element(size_t tag, size_t type, vector <Node> nodes)
{
    this->tag = tag;
    this->type = type;
    this->nodes = nodes;
}

// //

void make_shape()
{
    gmsh::model::add("t1_cube");

    double lc = 1e-2;

    // Make cube of points
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(0, 5, 0, lc, 2);
    gmsh::model::geo::addPoint(0, 0, 1, lc, 3);
    gmsh::model::geo::addPoint(0, 5, 1, lc, 4);

    gmsh::model::geo::addPoint(10, 0, 0, lc, 5);
    gmsh::model::geo::addPoint(10, 5, 0, lc, 6);
    gmsh::model::geo::addPoint(10, 0, 1, lc, 7);
    gmsh::model::geo::addPoint(10, 5, 1, lc, 8);

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

void generate_initial_mesh()
{
    double lc = 2; // 1 when actually testing

    // Add the first box (1st order corrections)
    gmsh::model::mesh::field::add("Box", 1);
    gmsh::model::mesh::field::setNumber(1, "VIn", lc / 8);
    gmsh::model::mesh::field::setNumber(1, "VOut", lc);
    gmsh::model::mesh::field::setNumber(1, "XMin", 1);
    gmsh::model::mesh::field::setNumber(1, "XMax", 10);
    gmsh::model::mesh::field::setNumber(1, "YMin", 1);
    gmsh::model::mesh::field::setNumber(1, "YMax", 4);
    gmsh::model::mesh::field::setNumber(1, "ZMin", 0);
    gmsh::model::mesh::field::setNumber(1, "ZMax", 1);

    // Add the second box (2nd order corrections)
    gmsh::model::mesh::field::add("Box", 2);
    gmsh::model::mesh::field::setNumber(2, "VIn", lc / 24);
    gmsh::model::mesh::field::setNumber(2, "VOut", lc);
    gmsh::model::mesh::field::setNumber(2, "XMin", 2);
    gmsh::model::mesh::field::setNumber(2, "XMax", 6);
    gmsh::model::mesh::field::setNumber(2, "YMin", 2);
    gmsh::model::mesh::field::setNumber(2, "YMax", 3);
    gmsh::model::mesh::field::setNumber(2, "ZMin", 0);
    gmsh::model::mesh::field::setNumber(2, "ZMax", 1);

    gmsh::model::mesh::field::add("Min", 3);
    gmsh::model::mesh::field::setNumbers(3, "FieldsList", {1,2});

    gmsh::model::mesh::field::setAsBackgroundMesh(3);


    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0); // 0 => mesh generation not based on Mesh.x
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    gmsh::model::mesh::generate(3);
}

vector<Node> get_mesh_nodes()
{
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords;
    std::vector<double> parametricCoords;

    gmsh::model::mesh::getNodesByElementType(4, nodeTags, coords, parametricCoords);

    vector<Node> nodes(nodeTags.size(), Node(0, vector<double>(3), vector<double>(3)));
    double x, y, z;
    double px, py, pz;

    for (unsigned i=0; i<nodeTags.size(); i++){
        x = coords[i*3];
        y = coords[i*3 + 1];
        z = coords[i*3 + 2];

//        cout << "parametricCoords.size() = " << parametricCoords.size() << '\n';
//        px = parametricCoords[i*3];
//        py = parametricCoords[i*3 + 1];
//        pz = parametricCoords[i*3 + 2];

        Node n(nodeTags[i], {x, y, z}, {px, py, pz});
        nodes[i] = n;
    }

    return nodes;
}

double target_size(double x, double y, double z)
{
//    return x / 10 + 0.1 ;
    if (x<1 && y<1) return 0.005;
    return 1;
//    return abs(cos(x)/2 * sin(y)/2) + 0.1;
}

//void remesh(int adapted_target_size_view_tag) {
//    gmsh::initialize();
//
//    make_shape();
//
//    try {
//        gmsh::merge("bg.pos");
//    } catch(...) {
//        gmsh::logger::write("Could not load background mesh: bye!", "error");
//        gmsh::finalize();
//    }
//
////    cout <<
//
//    make_shape();
//
//    // Add the post-processing view as a new size field:
//    int bg_field = gmsh::model::mesh::field::add("PostView");
//    gmsh::model::mesh::field::setNumber(bg_field, "ViewIndex", 0);
//
//    gmsh::model::mesh::field::setAsBackgroundMesh(bg_field);
//
//    // In order to compute the mesh sizes from the background mesh only, and
//    // disregard any other size constraints, one can set:
//    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0); // 0 => mesh generation not based on Mesh.x
//    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
//    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
//
//    gmsh::model::mesh::generate(3);
//
//    gmsh::write("reprocessed.msh");
//
//    gmsh::finalize();
//}

void get_elements(vector<Element> & elements)
{
    // Get elements
    std::vector<std::size_t> element_tags;
    std::vector<std::size_t> node_tags;
    gmsh::model::mesh::getElementsByType(4, element_tags, node_tags); // 4 = 4-node tetrahedron

    elements = vector<Element>(element_tags.size(), Element(0,0,{}));

    // Get Nodes
    vector<Node> nodes = get_mesh_nodes();

    // Evaluate target size
    for (unsigned i=0; i<nodes.size(); i++)
    {
        nodes[i].target_size = target_size(nodes[i].coord[0], nodes[i].coord[1], nodes[i].coord[2]);
    }

    // Make elements
//    vector<Element> elements(element_tags.size(), Element(0,0,{}));
    vector<Node> el_nodes(4, Node(0, {}, {}));

    for (unsigned i=0; i<elements.size(); i++)
    {
        for (unsigned j=0; j<4; j++)
        {
            el_nodes[j] = nodes[i*4+ j];
        }
        elements[i] = Element(element_tags[i], 4, el_nodes);
    }
}

// make background mesh
// Returns view tag
int add_refined_mesh_to_gmsh(vector<Element> &elements)
{
    // create a new view to be sued as the adapted target size mesh
    int view_tag = gmsh::view::add("adapted target size mesh");

//    vector<double> list_data(elements.size()*16);
    vector<double> list_data;

    for (unsigned e_i=0; e_i<elements.size(); e_i++)
    {
        Element e(elements[e_i]);

        vector<double> tet_coords(12); // 4 nodes * 3 coords (x1, ..., x4, y1, ...,y4, z1, ..., z4)
        vector<double> tet_target_sizes(4);
        for (unsigned i=0; i<e.nodes.size(); i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                tet_coords[i + 4*j] = e.nodes[i].coord[j];
            }
        }

        for (unsigned i=0; i<e.nodes.size(); i++)
        {
            tet_target_sizes[i] = e.nodes[i].target_size;
        }

        list_data.insert(list_data.end(), tet_coords.begin(), tet_coords.end());
        list_data.insert(list_data.end(), tet_target_sizes.begin(), tet_target_sizes.end());
    }

    gmsh::view::addListData(view_tag, "SS", elements.size(), list_data);

    int field_tag = gmsh::model::mesh::field::add("PostView");
    gmsh::model::mesh::field::setNumber(field_tag, "ViewTag", view_tag);

    gmsh::view::write(view_tag, "refined_bgmsh.pos");

    return field_tag;
}

vector<double> generate_initial_target_fields()
{
    double lc = 2;
    // Add the first box (1st order corrections)
    gmsh::model::mesh::field::add("Box", 1);
    gmsh::model::mesh::field::setNumber(1, "VIn", lc / 8);
    gmsh::model::mesh::field::setNumber(1, "VOut", lc);
    gmsh::model::mesh::field::setNumber(1, "XMin", 1);
    gmsh::model::mesh::field::setNumber(1, "XMax", 10);
    gmsh::model::mesh::field::setNumber(1, "YMin", 1);
    gmsh::model::mesh::field::setNumber(1, "YMax", 4);
    gmsh::model::mesh::field::setNumber(1, "ZMin", 0);
    gmsh::model::mesh::field::setNumber(1, "ZMax", 1);

    // Add the second box (2nd order corrections)
    gmsh::model::mesh::field::add("Box", 2);
    gmsh::model::mesh::field::setNumber(2, "VIn", lc / 24);
    gmsh::model::mesh::field::setNumber(2, "VOut", lc);
    gmsh::model::mesh::field::setNumber(2, "XMin", 2);
    gmsh::model::mesh::field::setNumber(2, "XMax", 6);
    gmsh::model::mesh::field::setNumber(2, "YMin", 2);
    gmsh::model::mesh::field::setNumber(2, "YMax", 3);
    gmsh::model::mesh::field::setNumber(2, "ZMin", 0);
    gmsh::model::mesh::field::setNumber(2, "ZMax", 1);

    return {(double)1, (double)2};
}

void mesh(vector<double> field_tags)
{
    int f_tag = gmsh::model::mesh::field::add("Min");
    gmsh::model::mesh::field::setNumbers(f_tag, "FieldsList", field_tags);

    gmsh::model::mesh::field::setAsBackgroundMesh(f_tag);

    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0); // 0 => mesh generation not based on Mesh.x
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    gmsh::model::mesh::generate(3);
}

int main()
{
    // initial mesh
    gmsh::initialize();

    make_shape();

    vector<double> bg_field_tags = generate_initial_target_fields();
    mesh(bg_field_tags);

    gmsh::write("basic_shape_initial.msh");

    vector<Element> elements;
    get_elements(elements);
    int number_initial_els = elements.size();

//    gmsh::clear()


    // refine mesh
    int refined_mesh_view_tag = add_refined_mesh_to_gmsh(elements);

    bg_field_tags.insert(bg_field_tags.end(), {(double)refined_mesh_view_tag});

    mesh(bg_field_tags);

    vector<Element> elements_ref;
    get_elements(elements_ref);

    int elements_ref_num = elements_ref.size();

    gmsh::write("refined.msh");

    gmsh::finalize();

    cout << "init: " << number_initial_els << ", refined: " << elements_ref_num << endl;
    cout << "refined_mesh_view_tag: " << refined_mesh_view_tag << endl;

    for (double x : bg_field_tags)
    {
        cout << x << ", ";
    }
    cout << endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    int num_el_nodes = 4; // 4 = 4 node tetrahedron
//
//    vector<Element> elements(element_tags.size(), Element(0, 0, {}));
//
//    vector<Node> el_nodes(num_el_nodes, Node(0, vector<double>(3), vector<double>(3)));
//    // For each node tag in each element, find the coordinates from vector<node> nodes, and calculate the target size for (x,y,z)
//    for (unsigned i=0; i<element_tags.size(); i++)
//    {
//        cout << i << endl;
//        for (unsigned j=0; j<num_el_nodes; j++)
//        {
//            el_nodes[j] = get_node(nodes, node_tags[i*num_el_nodes + j]);
//        }
//        elements[i] = Element(element_tags[i], 4, el_nodes);
//    }

    // create background mesh
    // remesh
//
//    gmsh::write("basic_shape_remesh.msh");

//    gmsh::finalize();

    return 0;
}

