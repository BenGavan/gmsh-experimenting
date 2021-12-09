//
// Created by Ben Gavan on 26/07/2021.
//

#include <iostream>
#include <fstream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;

double size_from_volume(double volume);

//const double constant = exp(-1.6383827341933301);
//const double N = 2.4659209630285632;

const double constant = exp(-1.3555);
const double N = 2.4245;

// Node //
class Node {
public:
    size_t tag;
    std::vector<double> coord;
    std::vector<double> parametric_coord;

    Node(std::size_t, std::vector<double>, std::vector<double>);
    void printCoords();
};

Node::Node(std::size_t t, std::vector<double> c, std::vector<double> pc)
{
    tag = t;
    coord = c;
    parametric_coord = pc;
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
//                double dv = levi_civita(i,j,k) * tet_vectors[i] * tet_vectors[3+j] * tet_vectors[6+k];
//                cout << "dv = " << dv << ", ";
//                cout << levi_civita(i,j,k) << ", ";
                volume += levi_civita(i,j,k) * tet_vectors[i] * tet_vectors[3+j] * tet_vectors[6+k];
            }
        }
    }
//    cout << endl;
    this->volume = abs(volume)/6.0;
//    cout << "Vol: " << volume << endl;
}

void Element::calculate_target_size_from_volume()
{
//    double constant = 1;
//    double N = 3;
    double l = pow(volume/constant, (double)1.0/N);
    linear_target_size = l;
}
// //


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

int add_uniform_target_size_field(double target_size=1)
{
    int tag = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(tag, "F", to_string(target_size));
    return tag;
}

void mesh_with_bg_fields(vector<double> field_tags)
{
    int t = gmsh::model::mesh::field::add("Min");
    gmsh::model::mesh::field::setNumbers(t, "FieldsList", field_tags);

    gmsh::model::mesh::field::setAsBackgroundMesh(t);

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
    vector<Node> el_nodes(4, Node(0, {}, {}));

    for (unsigned i=0; i<elements.size(); i++)
    {
        for (unsigned j=0; j<4; j++)
        {
            el_nodes[j] = nodes[i*4+ j];
        }
        elements[i] = Element(element_tags[i], 4, el_nodes);
    }

    for (unsigned i=0; i<elements.size(); i++)
    {
        for (unsigned j=0; j<4; j++)
        {
            elements[i].nodes[j].printCoords();
            cout << ", ";
        }
        cout << endl;
    }

    cout << "Nodes.size(): " << nodes.size() << endl;
    cout << "Elements.size(): " << elements.size() << endl;

    return elements;
}

void test_volume()
{
    Element e = Element(0, 0, {Node(0, {1+10, 0+10, 0+10}, {}),
                        Node(0, {1+10, 10+10, 10+10}, {}),
                        Node(0, {10+10, 1+10, 0+10}, {}),
                        Node(0, {0+10, 0+10, 10+10}, {})
    });

    e.calculate_volume();
    e.calculate_target_size_from_volume();
    cout << e.linear_target_size << endl;
}

void generate_average_size_data()
{


    double target_size;
    for (unsigned iteration=1; iteration<2; iteration++)
    {

        gmsh::initialize();
        gmsh::clear();

        make_cube();

//        target_size = .1 - 0.001*iteration;

        target_size = 1.0/(pow(iteration, 2));

        int uniform_field_tag = add_uniform_target_size_field(target_size);
        mesh_with_bg_fields({(double)uniform_field_tag});

        // - Find midpoint of elements
        // Get elements
        vector<Element> elements = get_elements();

        cout << "Number of elements: " << elements.size() << endl;

        ofstream vol_f;
        vol_f.open("volumes.txt", ios_base::app);

        // calculate volume of element tet
        for (unsigned i=0; i<elements.size(); i++)
        {
            // for info on tetrahedrons - https://people.sc.fsu.edu/~jburkardt/presentations/cg_lab_tetrahedrons.pdf
            elements[i].calculate_midpoint();
            elements[i].calculate_volume();
            elements[i].calculate_target_size_from_volume();
            elements[i].linear_target_size = size_from_volume(elements[i].volume);

            vol_f << target_size << "," << elements[i].linear_target_size << "," << elements[i].volume << '\n';
        }

        vol_f.close();

        cout << "--------------------\n--------------------\n--------------------\nTarget size: " << target_size << ", Number of elements: " << elements.size() << endl;

        double total_volume;
        double total_target_size = 0;
        for (unsigned i=0; i<elements.size(); i++)
        {
            total_volume += elements[i].volume;
            total_target_size += elements[i].linear_target_size;
        }
        double average_element_volume = total_volume/elements.size();
        double average_target_szie = total_target_size/elements.size();
        cout << "Average element volume: " << average_element_volume << endl;
        cout << "Average target size: " << average_target_szie << endl;

        ofstream f;
        f.open("out_.txt", ios_base::app);
//        f << target_size << "," << average_element_volume << '\n';
        f.close();

        gmsh::finalize();
    }



}

//const double constant = exp(-1.6383827341933301);
//const double N = 2.4659209630285632;

//const double constant = 1;
//const double N = 3;


double volume_from_size(double size)
{
    return constant * pow(size, N);
}

double size_from_volume(double volume)
{
    return pow(volume/constant, (double)1.0/N);
}

int main()
{
//    test_volume();
//    generate_average_size_data();
//    return 0;

    gmsh::initialize();

    // - Uniform mesh unit cube
    make_cube();

    int uniform_field_tag = add_uniform_target_size_field(0.0625);
    mesh_with_bg_fields({(double)uniform_field_tag});

    // - Find midpoint of elements
    // Get elements
    vector<Element> elements = get_elements();

    cout << "Number of elements: " << elements.size() << endl;

    // calculate volume of element tet
    double total_target_size = 0;
    double total_volume = 0;
    for (unsigned i=0; i<elements.size(); i++)
    {
        // for info on tetrahedrons - https://people.sc.fsu.edu/~jburkardt/presentations/cg_lab_tetrahedrons.pdf
        elements[i].calculate_midpoint();
        elements[i].calculate_volume();
        elements[i].calculate_target_size_from_volume();
//        cout << elements[i].linear_target_size << ", " << elements[i].volume << "\n";

        total_target_size += elements[i].linear_target_size;
        total_volume += elements[i].volume;
    }
    cout << endl;

    // calculate average target size
    double av_size = total_target_size/elements.size();
    double av_vol = total_volume/elements.size();

    cout << "Average size: " << av_size << endl;
    cout << "Average Volume: " << av_vol << endl;
    cout << "Total Volume: " << total_volume << endl;

    cout << "Size from average volume: " << size_from_volume(av_vol) << endl;
    cout << "Volume from average size: " << volume_from_size(av_size) << endl;


    // calculate what the target size is for that element volume (// TODO: Calculate target size (map the current element volume back to the target size))

    // - re-mesh using that size
    // make background mesh/field from elements


    // repeat and see if the size changes


    gmsh::finalize();

    return 0;
}