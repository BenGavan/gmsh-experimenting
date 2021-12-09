//
// Created by Ben Gavan on 08/08/2021.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include "../../libs/gmsh.h"

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

/// GridElement ///
class GridElement {
public:
    GridElement(int index, vector<double> min_corner, vector<double> max_corner) : index(index), min_corner(min_corner), max_corner(max_corner) {}

    vector<double> min_corner;
    vector<double> max_corner;
    vector<Element*> elements_pts;
    int index;

    double calculated_target_size = -1.;
    double target_size();
};

// naive target size of an element = average target size of mesh elements inside grid element
// (calculation of effective target size of grid element when there are no mesh elements inside grid element is
// delegated to Grid (since Grid has information of all grid elements) therefore
// Grid::target_size_for_grid_element(GridElement* grid_element_pt) should be used when using Grid as bg mesh).
double GridElement::target_size()
{
    if (elements_pts.size() == 0) return 0;
    double total_size = 0;
    for (unsigned i=0; i<elements_pts.size(); i++)
    {
        total_size += elements_pts[i]->linear_target_size;
    }
    return total_size/((double) elements_pts.size());
}

/// Grid //
class Grid {
public:
    Grid(vector<Element*> element_pts);

    vector<double> min_corner;
    vector<double> max_corner;
    vector<int> number_of_cells_for_axis = vector<int>(3); //  number of cells along a given axis
    vector<Element*> elements_pts;
    vector<GridElement*> grid_elements_pts; //  [(x_0,y_0,z_0),(x_1,y_0,z_0),...,(x_n,y_0,z_0),(x_0,y_1,z_0),(x_1,y_1,z_0),...,(x_n,y_n,z_0),(x_0,y_0,z_1),...,(x_n,y_n,z_n)]
    // define (x_0,y_0,z_0) as the origin and (x_0,y_0,z_0) = min_corner
    double average_mesh_element_size;

    void find_grid_corners();
    double average_element_volume();
    double average_element_size();
    double side_length(int side);
    double distance_from_origin_along_axis(int i, double p);
    void initialize_grid();
    int get_grid_element_index_for_coord(vector<double> coord);
    void process_elements();
    void calculate_grid_element_target_sizes();
//    double target_size_for_grid_element(GridElement* grid_element_pt);
    void calculate_average_element_size();
};

Grid::Grid(vector<Element*> element_pts) {
    this->elements_pts = element_pts;
    calculate_average_element_size();
    find_grid_corners();
    initialize_grid();
    process_elements();
}

void Grid::find_grid_corners()
{
    if (elements_pts.size() == 0) return;

    vector<double> mins(3);
    vector<double> maxs(3);

    for (unsigned i=0; i<3; i++)
    {
        mins[i] = elements_pts[0]->midpoint_coords[i];
        maxs[i] = elements_pts[0]->midpoint_coords[i];
    }

    for (unsigned i=1; i<elements_pts.size(); i++)
    {
        for (unsigned j=0; j<elements_pts[i]->nodes.size(); j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                if (elements_pts[i]->nodes[j].coord[k] < mins[k]) mins[k] = elements_pts[i]->nodes[j].coord[k];
                if (elements_pts[i]->nodes[j].coord[k] > maxs[k]) maxs[k] = elements_pts[i]->nodes[j].coord[k];
            }
        }
    }

    min_corner = mins;
    max_corner = maxs;

//    cout << "Max corner: " << to_string(max_corner) << endl;
//    cout << "Min corner: " << to_string(min_corner) << endl;
}

double Grid::average_element_volume()
{
    if (elements_pts.size()==0) return 0.;
    double total_volume = 0;
    for (unsigned i=0; i<elements_pts.size(); i++)
    {
        total_volume += elements_pts[i]->volume;
    }
    return total_volume / (double)elements_pts.size();
}

double Grid::average_element_size()
{
    if (elements_pts.size()==0) return 0.;
    double total_size = 0;
    for (unsigned i=0; i<elements_pts.size(); i++)
    {
        total_size += elements_pts[i]->linear_target_size;
    }
    return total_size / (double)elements_pts.size();
}

double Grid::side_length(int side)
{
    return max_corner[side] - min_corner[side];
}

void Grid::initialize_grid()
{
    double linear_size = average_mesh_element_size;
    //  number of elements along each axis
    int x_n = ceil(side_length(0) / linear_size);
    int y_n = ceil(side_length(1) / linear_size);
    int z_n = ceil(side_length(2) / linear_size);

    cout << "target size: " << linear_size << endl;
    cout << "average element size: " << average_element_size() << endl;

    cout << "x_n: " << x_n << endl;
    cout << "y_n: " << y_n << endl;
    cout << "z_n: " << z_n << endl;

    int num_element =  (x_n) * (y_n) * (z_n);
    this->grid_elements_pts = vector<GridElement*>(num_element);

    int y_counter = 0;
    int z_counter = 0;

    for (unsigned i=0; i<this->grid_elements_pts.size(); i++)
    {
        double min_x = (double)(i%x_n)*linear_size;
        double min_y = (double)(y_counter)*linear_size;
        double min_z = (double)z_counter*linear_size;

        double max_x = min_x + linear_size;
        double max_y = min_y + linear_size;
        double max_z = min_z + linear_size;

        vector<double> mins = {min_x, min_y, min_z};
        vector<double> maxs = {max_x, max_y, max_z};

//        cout << "Min: " << to_string(mins) << ", Max: " << to_string(maxs) << endl;

        if (i%x_n == x_n-1) y_counter++;
        if (y_counter == y_n)
        {
            y_counter = 0;
            z_counter++;
        }

        grid_elements_pts[i] = new GridElement(i, mins, maxs);
//        grid_elements_pts[i]->min_corner = mins;
//        grid_elements_pts[i]->max_corner = maxs;
//        grid_elements_pts[i]->index = i;
    }

    number_of_cells_for_axis = {x_n, y_n, z_n};
}

double Grid::distance_from_origin_along_axis(int i, double p)
{
    return p - min_corner[i];
}

int Grid::get_grid_element_index_for_coord(vector<double> coord)
{
    // total number of elements along each axis
    double target_size = this->average_mesh_element_size;

    int x_n = this->number_of_cells_for_axis[0];
    int y_n = this->number_of_cells_for_axis[1];
    int z_n = this->number_of_cells_for_axis[2];

    // number of elements along each axis the given coord should be in
//    int x_i = floor(distance_from_origin_along_axis(0, coord[0]) / target_size);
//    int y_i = floor(distance_from_origin_along_axis(1, coord[1]) / target_size);
//    int z_i = floor(distance_from_origin_along_axis(2, coord[2]) / target_size);

    int x_i = ceil(distance_from_origin_along_axis(0, coord[0]) / target_size) - 1;
    int y_i = ceil(distance_from_origin_along_axis(1, coord[1]) / target_size) - 1;
    int z_i = ceil(distance_from_origin_along_axis(2, coord[2]) / target_size) - 1;

    if (x_i < 0) x_i = 0;
    if (y_i < 0) y_i = 0;
    if (z_i < 0) z_i = 0;

    // [(x_0,y_0,z_0),(x_1,y_0,z_0),...,(x_n,y_0,z_0),(x_0,y_1,z_0),(x_1,y_1,z_0),...,(x_n,y_n,z_0),(x_0,y_0,z_1),...,(x_n,y_n,z_n)]

//    cout << target_size << endl;
    vector<double> v = {(double)x_i, (double) y_i, (double)z_i};
    int index = x_i + (y_i*x_n) + (z_i*x_n*y_n);
//    cout << "coord: " << to_string(coord) << ", pos: " << to_string(v) << ", index: " << index << endl;
//    print_vector(v);

    if (index > 1000000)
    {
        cout << "coord: " << to_string(coord) << ", pos: " << to_string(v) << ", index: " << index << endl;
        cout << "Target size: " << target_size << endl;
        cout << "dist 2:" << distance_from_origin_along_axis(2, coord[2]) << endl;
    }

    return x_i + (y_i*x_n) + (z_i*x_n*y_n);
}

void Grid::process_elements()
{
    initialize_grid();
    cout << "init grid then assign els" << endl;
    for (unsigned i=0; i<this->elements_pts.size(); i++)
    {
        int index = get_grid_element_index_for_coord(elements_pts[i]->midpoint_coords);
        this->grid_elements_pts[index]->elements_pts.push_back(elements_pts[i]);
    }
    cout << "Finished assigning mesh elements to grid" << endl;
}

// needs to be called before get_size_for_grid_element.
void Grid::calculate_grid_element_target_sizes()
{
    vector<GridElement*> zero_grid_elements_pts;

    // initial pass through all grid elements to calculate the target size if
    // the grid element contains one or more mesh elements
    for (unsigned i=0; i<this->grid_elements_pts.size(); i++)
    {
        double ts = grid_elements_pts[i]->target_size();
        if (ts > 0) {
            grid_elements_pts[i]->calculated_target_size = ts;
        } else {
            zero_grid_elements_pts.push_back(grid_elements_pts[i]);
        }
    }

//    cout << "Grid Min corner: " << to_string(this->min_corner) << ", Max: " << to_string(this->max_corner) << endl;

    // pass through the grid elements again
    // If grid element contains no mesh elements, take the average of the non-zero grid elements around it
    // repeat until all grid elements have been assigned a value
    while (zero_grid_elements_pts.size() > 0)
    {
        cout << "Number of zero elements: " << zero_grid_elements_pts.size() << endl;
        vector<GridElement*> remaining_zero_grid_elements_pts;
        for (unsigned i=0; i<zero_grid_elements_pts.size(); i++)
        {
            vector<int> indexes;

            // take average of surrounding elements
            int left_i = zero_grid_elements_pts[i]->index - 1;
            int right_i = zero_grid_elements_pts[i]->index + 1;

            int up_i = zero_grid_elements_pts[i]->index + this->number_of_cells_for_axis[0];
            int down_i = zero_grid_elements_pts[i]->index - this->number_of_cells_for_axis[0];

            int forward_i = zero_grid_elements_pts[i]->index + number_of_cells_for_axis[0]*number_of_cells_for_axis[1];
            int back_i = zero_grid_elements_pts[i]->index - number_of_cells_for_axis[0]*number_of_cells_for_axis[1];

//            if (i%number_of_cells_for_axis[0] != 0) indexes.push_back(left_i); // test if cell of on left edge (min of x)
//            if (i%number_of_cells_for_axis[0] != number_of_cells_for_axis[0]-1) indexes.push_back(right_i); // test is cell if on right edge (max x)

            if (zero_grid_elements_pts[i]->max_corner[0] < this->max_corner[0]) indexes.push_back(right_i); // test if cell if on right edge (max x)
            if (zero_grid_elements_pts[i]->min_corner[0] > this->min_corner[0]) indexes.push_back(left_i); // test if cell if on left edge (min x)

            if (zero_grid_elements_pts[i]->max_corner[1] < this->max_corner[1]) indexes.push_back(up_i); // test if cell if on top edge (max y)
            if (zero_grid_elements_pts[i]->min_corner[1] > this->min_corner[1]) indexes.push_back(down_i); // test if cell if on bottom edge (min y)

            if (zero_grid_elements_pts[i]->max_corner[2] < this->max_corner[2]) indexes.push_back(forward_i); // test if cell if on max front edge (depth) (max z)
            if (zero_grid_elements_pts[i]->max_corner[2] > this->max_corner[2]) indexes.push_back(back_i); // test if cell if on max front edge (depth) (max z)

            double total_size = 0;
            double non_zero_elements_included = 0;
            for (unsigned j=0; j<indexes.size(); j++)
            {
                if (grid_elements_pts[indexes[j]]->calculated_target_size <= 0) continue; // if this element is also zero, do not include in average
                total_size += this->grid_elements_pts[indexes[j]]->calculated_target_size;
                non_zero_elements_included++;
            }

//            cout << "indexes.size: " << indexes.size() << ", non zero elements includes: " << non_zero_elements_included << ", Total size: " << total_size << ", Min: " << to_string(zero_grid_elements_pts[i]->min_corner) << ", Max: " << to_string(zero_grid_elements_pts[i]->max_corner) << endl;

            if (non_zero_elements_included != 0) zero_grid_elements_pts[i]->calculated_target_size = total_size/non_zero_elements_included;
            if (zero_grid_elements_pts[i]->calculated_target_size <= 0) remaining_zero_grid_elements_pts.push_back(zero_grid_elements_pts[i]);
        }
        zero_grid_elements_pts = remaining_zero_grid_elements_pts;
    }
    cout << "Finished calculating target size for all gris elements (no more zero grid elements)" << endl;
}

void Grid::calculate_average_element_size()
{
    this->average_mesh_element_size = this->average_element_size();
}
/// END Grid ///


void make_cuboid(double x_length, double y_length, double z_length)
{
    gmsh::model::add("t1_cuboid");

    // Make cube of points
    int p1 = gmsh::model::geo::addPoint(0, 0, 0);
    int p2 = gmsh::model::geo::addPoint(0, y_length, 0);
    int p3 = gmsh::model::geo::addPoint(0, 0, z_length);
    int p4 = gmsh::model::geo::addPoint(0, y_length, z_length);

    int p5 = gmsh::model::geo::addPoint(x_length, 0, 0);
    int p6 = gmsh::model::geo::addPoint(x_length, y_length, 0);
    int p7 = gmsh::model::geo::addPoint(x_length, 0, z_length);
    int p8 = gmsh::model::geo::addPoint(x_length, y_length, z_length);

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

void make_cube(double side_length)
{
    make_cuboid(side_length, side_length, side_length);
}

int add_uniform_target_size_field(double target_size=1)
{
    int tag = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(tag, "F", to_string(target_size));
    return tag;
}

int add_background_box(vector<double> min, vector<double> max, double target_size)
{
    int t = gmsh::model::mesh::field::add("Box");
    gmsh::model::mesh::field::setNumber(t, "VIn", target_size);
    gmsh::model::mesh::field::setNumber(t, "VOut", target_size*100);
    gmsh::model::mesh::field::setNumber(t, "XMin", min[0]);
    gmsh::model::mesh::field::setNumber(t, "XMax", max[0]);
    gmsh::model::mesh::field::setNumber(t, "YMin", min[1]);
    gmsh::model::mesh::field::setNumber(t, "YMax", max[1]);
    gmsh::model::mesh::field::setNumber(t, "ZMin", min[2]);
    gmsh::model::mesh::field::setNumber(t, "ZMax", max[2]);
    return t;
}

vector<double> add_bg_fields_for_grid(Grid grid)
{
    grid.calculate_grid_element_target_sizes();

    uint64_t start_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    cout << "Adding background boxes: " << grid.grid_elements_pts.size() << endl;
    vector<double> tags;
    for (unsigned i=0; i<grid.grid_elements_pts.size(); i++)
    {
        double size = grid.grid_elements_pts[i]->calculated_target_size;
        double t = (double) add_background_box(grid.grid_elements_pts[i]->min_corner, grid.grid_elements_pts[i]->max_corner, size);
        tags.push_back(t);
        if (i%1000 == 0) {
            uint64_t end_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
            uint64_t duration = end_time - start_time;

            fstream f;
            string filename = "add_box_rate_0.07.txt";
            f.open(filename, ios_base::app);
            f << i << ',' << duration << '\n';
            f.close();

            cout << "Added box: " << i << ", Duration: " << duration / pow(10, 6) <<  endl;
        }
    }
    return tags;
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

    for (unsigned i=0; i<elements.size(); i++) {
        for (unsigned j = 0; j < 4; j++) {
            el_nodes[j] = nodes[i * 4 + j];
        }
        elements[i] = Element(element_tags[i], 4, el_nodes);
    }

    return elements;
}

vector<Element*> get_elements_pts()
{
    vector<Element> elements = get_elements();
    vector<Element*> element_pts(elements.size());
    for (unsigned i=0; i<elements.size(); i++)
    {
        element_pts[i] = &elements[i];
    }
    return element_pts;
}

vector<Element> make_uniform_mesh(double target_size)
{
    gmsh::initialize();
    gmsh::clear();

//    make_cube(1);
    make_cuboid(3,3,3);

    // Mesh
    int uniform_field_tag = add_uniform_target_size_field(target_size);
    mesh_with_bg_fields({(double)uniform_field_tag});

    vector<Element> elements = get_elements();

    gmsh::clear();
    gmsh::finalize();

    return elements;
}

void how_uniform_is_grid(Grid grid)
{
    for (unsigned i=0; i<grid.grid_elements_pts.size(); i++)
    {
        cout << "i: " << i << ", #elements: " << grid.grid_elements_pts[i]->elements_pts.size() << ", mid coords: ";
        for (unsigned j=0; j<grid.grid_elements_pts[i]->elements_pts.size(); j++)
        {
           cout << to_string(grid.grid_elements_pts[i]->elements_pts[j]->midpoint_coords) << ", ";
        }
        cout << endl;

        fstream f;
        f.open("grid_uniformity.txt", ios_base::app);
        f << i << ',' << grid.grid_elements_pts[i]->elements_pts.size() << '\n';
        f.close();
    }
}

void analyse_stability_iteration(int iteration, vector<Element*> element_pts, double target_size)
{
    double total_volume = 0;
    double total_target_size = 0;
    double elements_in_sample_volume = 0;

    // process mesh
    for (unsigned i=0; i<element_pts.size(); i++)
    {
        if (element_pts[i]->midpoint_coords[0] > .8 && element_pts[i]->midpoint_coords[0] < 2.2 &&
            element_pts[i]->midpoint_coords[1] > .8 && element_pts[i]->midpoint_coords[1] < 2.2 &&
            element_pts[i]->midpoint_coords[2] > .8 && element_pts[i]->midpoint_coords[2] < 2.2)
        {
            total_target_size += element_pts[i]->linear_target_size;
            total_volume += element_pts[i]->volume;
            elements_in_sample_volume += 1;
        }
    }

    double av_size = total_target_size/elements_in_sample_volume;
    double av_vol = total_volume/elements_in_sample_volume;

    cout << "Number of elements: " << element_pts.size() << endl;

    cout << "Average size: " << av_size << endl;
    cout << "Average Volume: " << av_vol << endl;
    cout << "Total Volume: " << total_volume << endl;

    cout << "Size from average volume: " << size_from_volume(av_vol) << endl;
    cout << "Volume from average size: " << volume_from_size(av_size) << endl;

    cout << "C = " << constant << ", N = " << N << endl;

    cout << "------------------------" << endl;

    string filename = "grid_stability_";
    filename.append(to_string(target_size));
    filename.append("_c=");
    filename.append(to_string(constant));
    filename.append("_N=");
    filename.append(to_string(N));
    filename.append(".txt");

    fstream f;
    f.open(filename, ios_base::app);
    f << iteration << ',' << av_size << ',' << av_vol << ',' << total_volume << ',' << elements_in_sample_volume << '\n';
    f.close();
}

int main()
{
    double target_size = 0.07;
    vector<Element> elements = make_uniform_mesh(target_size);
    for (unsigned i=0; i<elements.size(); i++)
    {
        elements[i].calculate_midpoint();
        elements[i].calculate_volume();
        elements[i].calculate_target_size_from_volume();
    }

    vector<Element*> element_pts;
    for (unsigned i=0; i<elements.size(); i++)
    {
        element_pts.push_back(&elements[i]);
    }

    cout << "Number of elements: " << elements.size() << endl;

    for (unsigned iter=0; iter<100; iter++)
    {
        // Analysis of previous iteration
        analyse_stability_iteration(iter, element_pts, target_size);

        cout << "**** iteration: " << iter << endl;

        Grid grid(element_pts);


        gmsh::initialize();
        gmsh::clear();

        // make shape
        make_cuboid(3, 3, 3);

        // make background mesh from grid
        vector<double> bg_fields = add_bg_fields_for_grid(grid);

        // mesh
        mesh_with_bg_fields(bg_fields);

        // extract elements
        elements = get_elements();
        element_pts = {};
        for (unsigned i=0; i<elements.size(); i++)
        {
            element_pts.push_back(&elements[i]);
        }

        for (unsigned i=0; i<element_pts.size(); i++)
        {
            element_pts[i]->calculate_midpoint();
            element_pts[i]->calculate_volume();
            element_pts[i]->calculate_target_size_from_volume();
        }

        gmsh::clear();
        gmsh::finalize();
    }



//    how_uniform_is_grid(grid);

    return 0;
}