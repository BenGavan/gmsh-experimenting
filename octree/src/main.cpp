//
// Created by Ben Gavan on 02/08/2021.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include "../libs/gmsh.h"

using namespace std;

const double constant = exp(-1.4891646142150403);
const double N = 2.9963431102916673;

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

//    this->calculate_midpoint();
//    this->calculate_volume();
//    this->calculate_target_size_from_volume();
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

    cout << "Nodes: " << nodes.size() << endl;

    for (unsigned i=0; i<nodes.size()-1; i++)
    {
        for (unsigned j=0; j<3; j++)
        {
            cout << "i: " << i << ", j: " << j << endl;
            cout << nodes.size() << endl;
            tet_vectors[i*3 + j] = nodes[i].coord[j] - nodes[3].coord[j];
        }
    }

    cout << to_string(tet_vectors) << endl;

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


void make_cuboid(double x_length, double y_length, double z_length)
{
    gmsh::model::add("t1_cuboid");

    double lc = 10;

    // Make cube of points
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(0, y_length, 0, lc, 2);
    gmsh::model::geo::addPoint(0, 0, z_length, lc, 3);
    gmsh::model::geo::addPoint(0, y_length, z_length, lc, 4);

    gmsh::model::geo::addPoint(x_length, 0, 0, lc, 5);
    gmsh::model::geo::addPoint(x_length, y_length, 0, lc, 6);
    gmsh::model::geo::addPoint(x_length, 0, z_length, lc, 7);
    gmsh::model::geo::addPoint(x_length, y_length, z_length, lc, 8);

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

void make_cube(double side_length)
{
    make_cuboid(side_length, side_length, side_length);
//    gmsh::model::add("t1_cube");
//
//    double lc = 10;
//
//    // Make cube of points
//    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
//    gmsh::model::geo::addPoint(0, side_length, 0, lc, 2);
//    gmsh::model::geo::addPoint(0, 0, side_length, lc, 3);
//    gmsh::model::geo::addPoint(0, side_length, side_length, lc, 4);
//
//    gmsh::model::geo::addPoint(side_length, 0, 0, lc, 5);
//    gmsh::model::geo::addPoint(side_length, side_length, 0, lc, 6);
//    gmsh::model::geo::addPoint(side_length, 0, side_length, lc, 7);
//    gmsh::model::geo::addPoint(side_length, side_length, side_length, lc, 8);
//
//    // Make lines to form cube
//    gmsh::model::geo::addLine(1, 2, 1);
//    gmsh::model::geo::addLine(1,3,2);
//    gmsh::model::geo::addLine(1,5,3);
//    gmsh::model::geo::addLine(2,6,4);
//    gmsh::model::geo::addLine(2,4,5);
//    gmsh::model::geo::addLine(6,5,6);
//    gmsh::model::geo::addLine(6,8,7);
//    gmsh::model::geo::addLine(8,4,8);
//    gmsh::model::geo::addLine(8,7,9);
//    gmsh::model::geo::addLine(4, 3,  10);
//    gmsh::model::geo::addLine(7, 3, 11);
//    gmsh::model::geo::addLine(5, 7, 12);
//
//    // Make curve loops to be able to make surfaces from them
//    gmsh::model::geo::addCurveLoop({1, 5, 10, -2}, 1);
//    gmsh::model::geo::addCurveLoop({4, 7, 8, -5}, 2);
//    gmsh::model::geo::addCurveLoop({1, 4, 6, -3}, 3);
//    gmsh::model::geo::addCurveLoop({3, 12, 11, -2}, 4);
//    gmsh::model::geo::addCurveLoop({6, 12, -9, -7}, 5);
//    gmsh::model::geo::addCurveLoop({8, 10, -11, -9}, 6);
//
//    // Make surfaces from curve loops
//    gmsh::model::geo::addPlaneSurface({1}, 1);
//    gmsh::model::geo::addPlaneSurface({2}, 2);
//    gmsh::model::geo::addPlaneSurface({3}, 3);
//    gmsh::model::geo::addPlaneSurface({4}, 4);
//    gmsh::model::geo::addPlaneSurface({5}, 5);
//    gmsh::model::geo::addPlaneSurface({6}, 6);
//
//    // Make surface loop from plane surfaces
//    gmsh::model::geo::addSurfaceLoop({1,2,3,4,5,6}, 1);
//
//    // Make volume from surface loop
//    gmsh::model::geo::addVolume({1}, 1);
//
//    gmsh::model::geo::synchronize();
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

vector<Element> make_uniform_mesh(double target_size)
{
    // Mesh
    int uniform_field_tag = add_uniform_target_size_field(target_size);
    mesh_with_bg_fields({(double)uniform_field_tag});

    vector<Element> elements = get_elements();

    return elements;
}

/// Octree ///

class OctreeNode
{
public:
    vector<double> min_corner;
    vector<double> max_corner;
    vector<Element*> elements;
    vector<OctreeNode*> child_node_pts;
    OctreeNode* parent_node_pt;
    int depth;
    bool is_refined = false;
    double target_linear_size;

    void refine();
    void calculate_corners();
    void calculate_base_node_corners();

    double mid(int i);
    int is_coord_in_upper(vector<double> c, int i);
    double side_length(int i);
    double average_volume();
    OctreeNode* get_base_node_pt();
    int get_max_depth();
};

double OctreeNode::mid(int i)
{
    return (max_corner[i] - min_corner[i]) / 2;
}

int OctreeNode::is_coord_in_upper(vector<double> c, int i)
{
    if (c[i] > mid(i)) return 1;
    return 0;
}

double OctreeNode::side_length(int i)
{
    return max_corner[i] - min_corner[i];
}

void OctreeNode::refine()
{
    this->is_refined = true;
    if (this->elements.size() < 4 || this->depth > this->get_base_node_pt()->get_max_depth())
    {
        this->is_refined = true;
        return;
    }

    // find the min/max corners of the 8 new sub octree nodes.
    if (child_node_pts.size() == 0)
    {
        for (unsigned i=0; i<8; i++)
        {
            OctreeNode* n = new OctreeNode();
            n->parent_node_pt = this;
            n->max_corner = vector<double>(3);
            n->min_corner = {};
            n->depth = this->depth + 1;

            double lower_x = this->min_corner[0] + (double)(1&i)*side_length(0)/2;
            double lower_y = this->min_corner[1] + (double)((2&i)/2)*side_length(1)/2;
            double lower_z = this->min_corner[2] + (double)((4&i)/4)*side_length(2)/2;

            n->min_corner = {lower_x, lower_y, lower_z};

            double upper_x = this->max_corner[0] + ((1&i)-1.)*side_length(0)/2;
            double upper_y = this->max_corner[1] + (((2&i)/2)-1.)*side_length(1)/2;
            double upper_z = this->max_corner[2] + (((4&i)/4)-1.)*side_length(2)/2;

            n->max_corner = {upper_x, upper_y, upper_z};

//            print_vector(n->min_corner);
//            print_vector(n->max_corner);
//            cout << endl;

//            for (unsigned j=0; j<3; j++)
//            {
//                n->min_corner[j] = this->min_corner[j] + (((1<<j)&i)/(1<<j))*side_length(j)/2;
//                n->max_corner[j] = this->max_corner[j] + (double)((((1<<j)&i)/(1<<j))-1.)*side_length(j)/2;
//            }

            child_node_pts.push_back(n);
        }
    }

    // assign elements to each of the new 8 child octree nodes.
    for (unsigned i=0; i<elements.size(); i++)
    {
        int pos = 4*is_coord_in_upper(elements[i]->midpoint_coords, 2)
                + 2*is_coord_in_upper(elements[i]->midpoint_coords, 1)
                + is_coord_in_upper(elements[i]->midpoint_coords, 0);

        child_node_pts[pos]->elements.push_back(elements[i]);
    }

    this->is_refined = true;
}

void OctreeNode::calculate_base_node_corners()
{
    if (elements.size() == 0) return;

    vector<double> mins(3);
    vector<double> maxs(3);

    for (unsigned i=0; i<3; i++)
    {
        mins[i] = elements[0]->midpoint_coords[i];
        maxs[i] = elements[0]->midpoint_coords[i];
    }

    for (unsigned i=1; i<elements.size(); i++)
    {
        for (unsigned j=0; j<elements[i]->nodes.size(); j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                if (elements[i]->nodes[j].coord[k] < mins[k]) mins[k] = elements[i]->nodes[j].coord[k];
                if (elements[i]->nodes[j].coord[k] > maxs[k]) maxs[k] = elements[i]->nodes[j].coord[k];
            }
        }
    }

    // find centre of mesh
    vector<double> centre(3);
    for (unsigned i=0; i<3; i++)
    {
        centre[i] = (mins[i] + maxs[i]) / 2;
    }

    // find max length
    double max_length = 0;
    for (unsigned i=0; i<3; i++)
    {
        double l = maxs[i] - mins[i];
        if (l > max_length) max_length = l;
    }

    cout << "Max dimensional length: " << max_length << endl;
    cout << "centre point: ";
    print_vector(centre);

    max_corner = vector<double>(3);
    min_corner = vector<double>(3);

    // cubify
    for (unsigned i=0; i<3; i++)
    {
        max_corner[i] = centre[i] + max_length/2;
        min_corner[i] = centre[i] - max_length/2;
    }
    cout << "Max corner: " << to_string(max_corner) << endl;
    cout << "Min corner: " << to_string(min_corner) << endl;

//    min_corner = mins;
//    max_corner = maxs;
}

double OctreeNode::average_volume()
{
    double total_volume = 0;
    for (unsigned i=0; i<elements.size(); i++)
    {
        total_volume += elements[i]->volume;
    }
    return total_volume / (double)elements.size();
}

OctreeNode* OctreeNode::get_base_node_pt()
{
    if (this->depth == 0) return this;
    return this->parent_node_pt->get_base_node_pt();
}

int OctreeNode::get_max_depth()
{
    double num = log(this->get_base_node_pt()->elements.size());
    double den = log(8);
    return num/den;
}

/// END Octree Class ///

//vector<OctreeNode*> refine_recursive(vector<OctreeNode*> o_node_pts)
//{
//    vector<OctreeNode*> child_node_pts;
//    for (unsigned i=0; i<o_node_pts.size(); i++)
//    {
//        o_node_pts[i]->refine();
//        if (!o_node_pts[i]->child_node_pts.empty())
//        {
//            child_node_pts.insert(child_node_pts.end(), o_node_pts[i]->child_node_pts.begin(), o_node_pts[i]->child_node_pts.end());
//        }
//    }
//    if (child_node_pts.empty()) return o_node_pts;
//    return refine_recursive(child_node_pts);
//}

// returns the new refined nodes as pointers.
vector<OctreeNode*> refine_nodes(vector<OctreeNode*> o_n_pts)
{
    vector<OctreeNode*> new_node_pts;
    for (unsigned i=0; i<o_n_pts.size(); i++)
    {
        if (!o_n_pts[i]->is_refined && o_n_pts[i]->child_node_pts.size() == 0)
        {
            o_n_pts[i]->refine();
            new_node_pts.insert(new_node_pts.end(), o_n_pts[i]->child_node_pts.begin(), o_n_pts[i]->child_node_pts.end());
        }
    }
    return new_node_pts;
}

//vector<OctreeNode*> get_leaves(vector<OctreeNode*> oct_nodes_pts)
//{
//    vector<OctreeNode*> leaves_pts;
//    for (unsigned i=0; i<oct_nodes_pts.size(); i++)
//    {
//        if (oct_nodes_pts.)
//        leaves_pts.push_back()
//    }
//}

int add_background_box(vector<double> min, vector<double> max, double target_size)
{
    int t = gmsh::model::mesh::field::add("Box");
    gmsh::model::mesh::field::setNumber(1, "VIn", target_size);
    gmsh::model::mesh::field::setNumber(1, "VOut", 1);
    gmsh::model::mesh::field::setNumber(1, "XMin", min[0]);
    gmsh::model::mesh::field::setNumber(1, "XMax", max[0]);
    gmsh::model::mesh::field::setNumber(1, "YMin", min[1]);
    gmsh::model::mesh::field::setNumber(1, "YMax", max[1]);
    gmsh::model::mesh::field::setNumber(1, "ZMin", min[2]);
    gmsh::model::mesh::field::setNumber(1, "ZMax", max[2]);
    return t;
}

// returns: field tags (vector<double>) (int->double handled)
vector<double> octree_to_background_mesh(vector<OctreeNode*> oct_nodes_pts)
{
    vector<double> tags(oct_nodes_pts.size());
    for (unsigned i=0; i<oct_nodes_pts.size(); i++)
    {
        double target_size = size_from_volume(oct_nodes_pts[i]->average_volume());
        oct_nodes_pts[i]->target_linear_size = target_size;

        tags[i] = add_background_box(oct_nodes_pts[i]->min_corner, oct_nodes_pts[i]->max_corner, target_size);
    }
    return tags;
}

vector<OctreeNode*> tree_from_els(vector<Element*> element_pts)
{
    OctreeNode base_node;
    base_node.elements = element_pts;
    base_node.calculate_base_node_corners();
    base_node.depth = 0;

    // make tree by node refinement
    vector<OctreeNode*> all_nodes_pts = {&base_node};
    vector<OctreeNode*> new_nodes_pts;

    int x = 10;
    do {
        new_nodes_pts = refine_nodes(all_nodes_pts);
        all_nodes_pts.insert(all_nodes_pts.end(), new_nodes_pts.begin(), new_nodes_pts.end());
        x--;

        for (unsigned i=0; i<all_nodes_pts.size(); i++)
        {
//            cout << all_nodes_pts[i] << " " << all_nodes_pts[i]->is_refined << " " << all_nodes_pts[i]->elements.size() << " " << all_nodes_pts[i]->depth << " " << all_nodes_pts[i]->parent_node_pt;
//            for (unsigned j=0; j<all_nodes_pts[i]->elements.size(); j++)
//            {
//                cout << to_string(all_nodes_pts[i]->elements[j]->midpoint_coords);
//            }
//            cout << " Min: " << to_string(all_nodes_pts[i]->min_corner) << ", Max: " << to_string(all_nodes_pts[i]->max_corner) << " " << all_nodes_pts[i]->max_corner[0] << endl;

//            cout << " " << all_nodes_pts.size() << endl;
//            cout << all_nodes_pts[i]->parent_node_pt << endl;
            cout << all_nodes_pts[i]->get_base_node_pt()  << ", #elements: " << all_nodes_pts[i]->elements.size() << endl;
        }
        cout << "----- " << x << endl;
//        for (unsigned i=0; i<new_nodes_pts.size(); i++)
//        {
//            cout << to_string(new_nodes_pts[i]->min_corner) << ", ";
//        }
//        cout << endl;
    } while (new_nodes_pts.size() != 0);

    return all_nodes_pts;
}

vector<Element*> elements_to_elements_pts(vector<Element> &elements)
{
    vector<Element*> elements_pts = vector<Element*>(elements.size());
    for (unsigned i=0; i<elements.size(); i++)
    {
        elements_pts[i] = &elements[i];
    }
    return elements_pts;
}

void append_to_file(int i, int number_elements, double average_size, double average_volume, double total_volume)
{
    fstream f;
    f.open("out_0.03.txt", ios_base::app);
    f << i << ',' << number_elements << ',' << average_size << ',' << average_volume << ',' << total_volume << '\n';
    f.close();
}

void analyse(vector<Element*> elements_pts, int iteration)
{
    double total_target_size = 0;
    double total_volume = 0;
    double elements_in_test_volume = 0;
    for (unsigned i=0; i<elements_pts.size(); i++)
    {
        elements_pts[i]->calculate_midpoint();
        elements_pts[i]->calculate_volume();
        elements_pts[i]->calculate_target_size_from_volume();

        if (elements_pts[i]->midpoint_coords[0] > .7 && elements_pts[i]->midpoint_coords[0] < 2.3 &&
            elements_pts[i]->midpoint_coords[1] > .7 && elements_pts[i]->midpoint_coords[1] < 2.3 &&
            elements_pts[i]->midpoint_coords[2] > .7 && elements_pts[i]->midpoint_coords[2] < 2.3)
        {
            total_target_size += elements_pts[i]->linear_target_size;
            total_volume += elements_pts[i]->volume;
            elements_in_test_volume += 1;
        }
    }

    // calculate average target size
    double av_size = total_target_size/elements_in_test_volume;
    double av_vol = total_volume/elements_in_test_volume;

    append_to_file(iteration, elements_in_test_volume, av_size, av_vol, total_volume);
}

void process_elements_pts(vector<Element*> elements_pts)
{
    for (unsigned i=0; i<elements_pts.size(); i++)
    {
        elements_pts[i]->calculate_midpoint();
        elements_pts[i]->calculate_volume();
        elements_pts[i]->calculate_target_size_from_volume();
    }
}

int main()
{
    gmsh::initialize();
    gmsh::clear();

    make_cube(3);
//    make_cuboid(3,3,3);

    vector<Element> elements = make_uniform_mesh(1);

    gmsh::clear();
    gmsh::finalize();

//    vector<Element*> elements_pts = elements_to_elements_pts(elements);

    vector<Element*> elements_pts;
    for (unsigned i=0; i<elements.size(); i++)
    {
        elements_pts.push_back(&elements[i]);
    }

    process_elements_pts(elements_pts);

    cout << "Seg fault: " << endl;
    analyse(elements_pts, -1);
    cout << "Not here" << endl;

    cout << "Number of elements: " << elements_pts.size() << endl;

    vector<OctreeNode*> tree = tree_from_els(elements_pts);

    // stability
    for (unsigned iter=0; iter<100; iter++)
    {

        gmsh::initialize();
        gmsh::clear();

        make_cuboid(3,3,3);
        vector<double> bg_tags = octree_to_background_mesh(tree);
        mesh_with_bg_fields(bg_tags);

        vector<Element> elements = get_elements();

        elements_pts = elements_to_elements_pts(elements);
        cout << "Seg fault: " << endl;
        analyse(elements_pts, iter);
        cout << "Not here" << endl;

        tree = tree_from_els(elements_pts);

        gmsh::clear();
        gmsh::finalize();
    }

    return 0;
}
