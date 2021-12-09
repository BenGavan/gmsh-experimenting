//
// Created by Ben Gavan on 30/07/2021.
// for info on tetrahedrons - https://people.sc.fsu.edu/~jburkardt/presentations/cg_lab_tetrahedrons.pdf
//

#include <iostream>
#include <fstream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;

//const double constant = exp(-1.4891646142150403); // = 0.22556100731
//const double constant = 0.23210085272629383;
const double constant = 0.2215;
//const double N = 2.9963431102916673;
const double N = 3.0;

double volume_from_size(double size)
{
    return constant * pow(size, N);
}

double size_from_volume(double volume)
{
    return pow(volume/constant, (double)1.0/N);
}


// Node //
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
// //


void make_cube(double side_length)
{
    gmsh::model::add("t1_cube");

    double lc = 10;

    // Make cube of points
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(0, side_length, 0, lc, 2);
    gmsh::model::geo::addPoint(0, 0, side_length, lc, 3);
    gmsh::model::geo::addPoint(0, side_length, side_length, lc, 4);

    gmsh::model::geo::addPoint(side_length, 0, 0, lc, 5);
    gmsh::model::geo::addPoint(side_length, side_length, 0, lc, 6);
    gmsh::model::geo::addPoint(side_length, 0, side_length, lc, 7);
    gmsh::model::geo::addPoint(side_length, side_length, side_length, lc, 8);

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

    for (unsigned i=0; i<elements.size(); i++)
    {
        for (unsigned j=0; j<4; j++)
        {
            el_nodes[j] = nodes[i*4+ j];
        }
        elements[i] = Element(element_tags[i], 4, el_nodes);
    }

    return elements;
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
            tet_target_sizes[i] = e.linear_target_size;
        }

        list_data.insert(list_data.end(), tet_coords.begin(), tet_coords.end());
        list_data.insert(list_data.end(), tet_target_sizes.begin(), tet_target_sizes.end());
    }

    gmsh::view::addListData(view_tag, "SS", elements.size(), list_data);

    int field_tag = gmsh::model::mesh::field::add("PostView");
    gmsh::model::mesh::field::setNumber(field_tag, "ViewTag", view_tag);

//    gmsh::view::write(view_tag, "refined_bgmsh.pos");

    return field_tag;
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

// create initial uniform mesh
vector<Element> initial_uniform_mesh(double target_size)
{
    gmsh::initialize();
    gmsh::clear();

    make_cube(3);

    // uniform mesh
    int uniform_field_tag = add_uniform_target_size_field(target_size);
    mesh_with_bg_fields({(double)uniform_field_tag});

    // extract mesh
    vector<Element> elements = get_elements();

    gmsh::clear();
    gmsh::finalize();

    return elements;
}

int main()
{
    cout << "C = " << constant << ", N = " << N << endl;

    uint64_t start_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // create initial uniform mesh
    vector<Element> elements = initial_uniform_mesh(0.03);

    for (unsigned iter=0; iter<100; iter++)
    {

        double total_volume = 0;
        double total_target_size = 0;
        double elements_in_sample_volume = 0;

        // process mesh
        for (unsigned i=0; i<elements.size(); i++)
        {
            elements[i].calculate_midpoint();
            elements[i].calculate_volume();
            elements[i].calculate_target_size_from_volume();

            if (elements[i].midpoint_coords[0] > 1 && elements[i].midpoint_coords[0] < 2 &&
                elements[i].midpoint_coords[1] > 1 && elements[i].midpoint_coords[1] < 2 &&
                elements[i].midpoint_coords[2] > 1 && elements[i].midpoint_coords[2] < 2)
            {
                total_target_size += elements[i].linear_target_size;
                total_volume += elements[i].volume;
                elements_in_sample_volume += 1;
            }
        }

        double av_size = total_target_size/elements_in_sample_volume;
        double av_vol = total_volume/elements_in_sample_volume;uint64_t end_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        uint64_t duration = end_time - start_time;



        cout << "Number of elements: " << elements.size() << endl;

        cout << "Average size: " << av_size << endl;
        cout << "Average Volume: " << av_vol << endl;
        cout << "Total Volume: " << total_volume << endl;

        cout << "Size from average volume: " << size_from_volume(av_vol) << endl;
        cout << "Volume from average size: " << volume_from_size(av_size) << endl;

        cout << "C = " << constant << ", N = " << N << endl;

        cout << "Duration: " << duration << endl;

        cout << "------------------------" << endl;

        fstream f;
        f.open("out_cheat_mesh_stability_0.03_c=0.2215_N=3.txt", ios_base::app);
        f << iter << ',' << av_size << ',' << av_vol << ',' << total_volume << ',' << elements_in_sample_volume << ',' << duration << '\n';
        f.close();

        // new mesh

        gmsh::initialize();
        gmsh::clear();

        // make new target background field from elements[i].linear_target_size

        make_cube(3);
        int refined_filed_tag = add_refined_mesh_to_gmsh(elements);

        mesh({(double)refined_filed_tag});

        elements = get_elements();

        gmsh::clear();
        gmsh::finalize();

        start_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    }

    return 0;
}