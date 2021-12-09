//
// Created by Ben Gavan on 11/08/2021.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include "../libs/gmsh.h"

using namespace std;

//const double constant = exp(-1.4891646142150403);
const double constant = 0.222;
//const double N = 2.9963431102916673;
const double N = 3;

double volume_from_size(double size)
{
    return constant * pow(size, N);
}

double size_from_volume(double volume)
{
    return pow(volume/constant, (double)1.0/N);
}

string to_string(vector<double> v)
{
    string s;
    s += '(';
    for (unsigned i=0; i<v.size(); i++)
    {
        s += to_string(v[i]);
        if (i != v.size()-1) s += ", ";
    }
    s += ")";
    return s;
}

void print_vector(vector<double> v)
{
    cout << to_string(v) << endl;
}


/// Node ///
class Node {
public:
    size_t tag;
    std::vector<double> coord;

    Node(std::size_t, std::vector<double>);
    void printCoords();
};

Node::Node(std::size_t t, std::vector<double> c)
{
    tag = t;
    coord = c;
}

void Node::printCoords()
{
    cout << "(" << coord[0] << ", " << coord[1] << ", " << coord[2] << ")";
}
/// ///

// Element //
class Element
{
public:
    size_t tag;
    size_t type;
    vector<Node> nodes;
    vector<double> midpoint_coords;
    double volume;
    double linear_target_size;

    Element(size_t tag, size_t type, vector<Node> nodes);
    void calculate_midpoint();
    void calculate_volume();
    void calculate_target_size_from_volume();
};

Element::Element(size_t tag, size_t type, vector <Node> nodes)
{
    this->tag = tag;
    this->type = type;
    this->nodes = nodes;
}

void Element::calculate_midpoint()
{
    // midpoint = centre of mass
    vector<double> com_coords(3);
    for (unsigned i=0; i<nodes.size(); i++)
    {
        for (unsigned j=0; j<3; j++)
        {
            com_coords[j] += nodes[i].coord[j] / nodes.size();
        }
    }
    midpoint_coords = com_coords;
}

// = 1 for i!=j!=k and
int levi_civita(int i, int j, int k)
{
    return 0.5 * (i-j) * (j-k) * (k-i);
}

void Element::calculate_volume()
{
    // Volume of tetrahedron = (1/6)|axb.c| = (1/6)|v_1 x v_2 . v_3|
    // find the defining vectors a,b,c from the 4 global positions from the node global position vectors
    // define node 4 (p_4) to be the origin of the tetrahedron.
    // v_i = p_i - p_4

    vector<double> tet_vectors(3*3); // [a_1,...,a_3, b_1,...b_3, c_1,...c_3]

    for (unsigned i=0; i<nodes.size()-1; i++)
    {
        for (unsigned j=0; j<3; j++)
        {
            tet_vectors[i*3 + j] = nodes[i].coord[j] - nodes[3].coord[j];
        }
    }
//    cout << "(";
//    for (unsigned i=0; i<tet_vectors.size(); i++)
//    {
//        cout << tet_vectors[i];
//        if (i != tet_vectors.size()-1) cout << ", ";
//    }
//    cout << ")" << endl;

    // calculate volume
    this->volume = 0;
    for (unsigned i=0; i<3; i++)
    {
        for (unsigned j=0; j<3; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                volume +=  levi_civita(i,j,k) * tet_vectors[i] * tet_vectors[3+j] * tet_vectors[6+k];
            }
        }
    }
    this->volume = abs(volume)/6.0;
}

void Element::calculate_target_size_from_volume()
{
    linear_target_size = pow(volume/constant, (double)1.0/N);
}
/// END Element ///

vector<Node> get_mesh_nodes()
{
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords;
    std::vector<double> parametricCoords;

    gmsh::model::mesh::getNodesByElementType(4, nodeTags, coords, parametricCoords);

    vector<Node> nodes(nodeTags.size(), Node(0, vector<double>(3)));
    double x, y, z;

    for (unsigned i=0; i<nodeTags.size(); i++){
        x = coords[i*3];
        y = coords[i*3 + 1];
        z = coords[i*3 + 2];

        Node n(nodeTags[i], {x, y, z});
        nodes[i] = n;
    }

    return nodes;
}

vector<Element> get_elements()
{
    // Get elements
    vector<Element> elements;
    std::vector<std::size_t> element_tags;
    std::vector<std::size_t> node_tags;
    gmsh::model::mesh::getElementsByType(4, element_tags, node_tags); // 4 = 4-node tetrahedron

    elements = vector<Element>(element_tags.size(), Element(0,0,{}));

    // Get Nodes
    vector<Node> nodes = get_mesh_nodes();

    // Make elements
    vector<Node> el_nodes(4, Node(0, {}));

    for (unsigned i=0; i<elements.size(); i++) {
        for (unsigned j = 0; j < 4; j++) {
            el_nodes[j] = nodes[i * 4 + j];
        }
        elements[i] = Element(element_tags[i], 4, el_nodes);
    }

    return elements;
}

void make_cuboid(vector<double> origin_coords, vector<double> lengths)
{

}

void make_cube()
{
    gmsh::model::add("cube");

    // Make cube of points
    int p1 = gmsh::model::geo::addPoint(0, 0, 0);
    int p2 = gmsh::model::geo::addPoint(0, 1, 0);
    int p3 = gmsh::model::geo::addPoint(0, 0, 1);
    int p4 = gmsh::model::geo::addPoint(0, 1, 1);

    int p5 = gmsh::model::geo::addPoint(1, 0, 0);
    int p6 = gmsh::model::geo::addPoint(1, 1, 0);
    int p7 = gmsh::model::geo::addPoint(1, 0, 1);
    int p8 = gmsh::model::geo::addPoint(1, 1, 1);

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

    gmsh::model::geo::synchronize();
}

int add_background_box(vector<double> min, vector<double> max, double target_size)
{
    int t = gmsh::model::mesh::field::add("Box");
    gmsh::model::mesh::field::setNumber(t, "VIn", target_size);
    gmsh::model::mesh::field::setNumber(t, "VOut", 1);
    gmsh::model::mesh::field::setNumber(t, "XMin", min[0]);
    gmsh::model::mesh::field::setNumber(t, "XMax", max[0]);
    gmsh::model::mesh::field::setNumber(t, "YMin", min[1]);
    gmsh::model::mesh::field::setNumber(t, "YMax", max[1]);
    gmsh::model::mesh::field::setNumber(t, "ZMin", min[2]);
    gmsh::model::mesh::field::setNumber(t, "ZMax", max[2]);
    return t;
}

vector<double> add_bg_fields()
{
    vector<double> tags;

    double size_1 = 0.1;
    vector<double> min_1 = {0.1,0.1,0.1};
//    vector<double> max_1 = {1.717644, 1.717644, 1.717644};
    vector<double> max_1 = {.2, .2, .2};

    double t1 = (double) add_background_box(min_1, max_1, size_1);

    tags.push_back(t1);

    double t2 = (double) add_background_box({.15, .15, .15}, {.17, .17, .17}, 0.05);

    tags.push_back(t2);

    int t_e = gmsh::model::mesh::field::add("Box");

    tags.push_back(t_e);

    print_vector(tags);

    return tags;
}

void mesh_with_bg_fields(vector<double> field_tags)
{

    int f_tag = gmsh::model::mesh::field::add("Min");
    gmsh::model::mesh::field::setNumbers(f_tag, "FieldsList", field_tags);

    gmsh::model::mesh::field::setAsBackgroundMesh(f_tag);

    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0); // 0 => mesh generation not based on Mesh.x
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    gmsh::model::mesh::generate(3);
//    int t = gmsh::model::mesh::field::add("Min");
//    gmsh::model::mesh::field::setNumbers(t, "FieldsList", field_tags);
//
//    gmsh::model::mesh::field::setAsBackgroundMesh(t);
//
//    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0); // 0 => mesh generation not based on Mesh.x
//    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
//    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
//
//    gmsh::model::mesh::generate(3);
}



int main()
{
    gmsh::initialize();
    gmsh::clear();

    // make shape
    make_cube();

    // mesh
    vector<double> tags = add_bg_fields();
    mesh_with_bg_fields(tags);

    vector<Element> elements = get_elements();
    cout << "Number of elements: " << elements.size() << endl;


    gmsh::clear();
    gmsh::finalize();

    return 0;
}

