//
// Created by Ben Gavan on 08/08/2021.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include "../libs/gmsh.h"

using namespace std;

//const double constant = exp(-1.4891646142150403); // from (naive) best fit
//const double constant = 0.2152; // this is what's currently best
double constant = 0.2152; // not const since it is the variable to be are adjusting to find the value for each target size.
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

string to_string(vector<int> v)
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

void print_vector(vector<int> v)
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
    vector<double> cell_size = vector<double>(3); // size of each cubic cell which makes up the grid (dx, dy, dz)
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
    void write(string filename);
    Grid shrink_by(double factor); // shrinks the elements being analysed inwards by a factor (x_min+factor, x_max-factor), returns a new grid
};

Grid Grid::shrink_by(double factor)
{
    vector<Element*> new_grid_elements;
    vector<double> new_min_corner(3);
    vector<double> new_max_corner(3);

    for (unsigned i=0; i<3; i++)
    {
        new_min_corner[i] = min_corner[i] + factor;
        new_max_corner[i] = max_corner[i] - factor;
    }

    for (unsigned i=0; i<elements_pts.size(); i++)
    {
        if ((elements_pts[i]->midpoint_coords[0] > new_min_corner[0] || elements_pts[i]->midpoint_coords[0] < new_max_corner[0])
         && (elements_pts[i]->midpoint_coords[1] > new_min_corner[1] || elements_pts[i]->midpoint_coords[1] < new_max_corner[1])
         && (elements_pts[i]->midpoint_coords[2] > new_min_corner[2] || elements_pts[i]->midpoint_coords[2] < new_max_corner[2]))
        {
            new_grid_elements.push_back(elements_pts[i]);
        }
    }

    return Grid(new_grid_elements);
}

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

    double dx = side_length(0) / double(x_n);
    double dy = side_length(1) / double(y_n);
    double dz = side_length(2) / double(z_n);

    this->cell_size = {dx, dy, dz};
    this->number_of_cells_for_axis = {x_n, y_n, z_n};

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
        double min_x = this->min_corner[0] + (double)(i%x_n)*dx;
        double min_y = this->min_corner[1] + (double)(y_counter)*dy;
        double min_z = this->min_corner[2] + (double)z_counter*dz;

        double max_x = min_x + dx;
        double max_y = min_y + dy;
        double max_z = min_z + dz;

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
    }

    number_of_cells_for_axis = {x_n, y_n, z_n};
}

double Grid::distance_from_origin_along_axis(int i, double p)
{
    return p - min_corner[i];
}

int Grid::get_grid_element_index_for_coord(vector<double> coord)
{
    // Structure of grid element vector:
    // [(x_0,y_0,z_0),(x_1,y_0,z_0),...,(x_n,y_0,z_0),(x_0,y_1,z_0),(x_1,y_1,z_0),...,(x_n,y_n,z_0),(x_0,y_0,z_1),...,(x_n,y_n,z_n)]

    // total number of elements along each axis
    double dx = this->cell_size[0];
    double dy = this->cell_size[1];
    double dz = this->cell_size[2];

    int x_n = this->number_of_cells_for_axis[0];
    int y_n = this->number_of_cells_for_axis[1];
    int z_n = this->number_of_cells_for_axis[2];

    // number of elements along each axis the given coord should be in
//    int x_i = floor(distance_from_origin_along_axis(0, coord[0]) / target_size);
//    int y_i = floor(distance_from_origin_along_axis(1, coord[1]) / target_size);
//    int z_i = floor(distance_from_origin_along_axis(2, coord[2]) / target_size);

    int x_i = ceil(distance_from_origin_along_axis(0, coord[0]) / dx) - 1;
    int y_i = ceil(distance_from_origin_along_axis(1, coord[1]) / dy) - 1;
    int z_i = ceil(distance_from_origin_along_axis(2, coord[2]) / dz) - 1;

    if (x_i < 0) x_i = 0;
    if (y_i < 0) y_i = 0;
    if (z_i < 0) z_i = 0;

    vector<double> v = {(double)x_i, (double) y_i, (double)z_i};
    int index = x_i + (y_i*x_n) + (z_i*x_n*y_n);

    return x_i + (y_i*x_n) + (z_i*x_n*y_n);
}

void Grid::process_elements()
{
    initialize_grid();
    cout << "Initialize grid then assign els" << endl;
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
            // for now just assign the average target size to zero grid element elements
            grid_elements_pts[i]->calculated_target_size = average_mesh_element_size;
//            zero_grid_elements_pts.push_back(grid_elements_pts[i]);
        }
    }

    return;

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

int add_step_uniform_target_size_field()
{
    int tag = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(tag, "F", "0.02 + 0.003*x");
    return tag;
}

int add_background_box(vector<double> min, vector<double> max, double target_size)
{
    int t = gmsh::model::mesh::field::add("Box");
    gmsh::model::mesh::field::setNumber(t, "VIn", target_size);
    gmsh::model::mesh::field::setNumber(t, "VOut", 0.1);
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

    // create a new view to be sued as the adapted target size mesh
    int view_tag = gmsh::view::add("adapted target size mesh");

    cout << "Adding background Mesh with grid elements: " << grid.grid_elements_pts.size() << endl;

    vector<double> list_data;

    //    = {0,1,0,0,1,0, 0,0,0,1,1,1, 0,1,1,0,1,1, .5,.5,.5,.5,.5,.5,
//                                0,1,1,0,1,1, 0,0,0,1,1,1, 0,1,0,0,1,0, 5,5,5,5,5,5};
//    vector<double> l2 = {0,1,1,0,1,1, 0,0,0,1,1,1, 0,1,0,0,1,0, .5,.5,.5,.5,.5,.5};

//        SI(0,0,0, 1,0,1, 1,0,0, 0,1,0, 1,1,1, 1,1,0){.5,.5,.5,.5,.5,.5};

    for (unsigned i=0; i<grid.grid_elements_pts.size(); i++)
    {
        double size = grid.grid_elements_pts[i]->calculated_target_size;

        vector<double> p1_coords = {
                grid.grid_elements_pts[i]->min_corner[0],
                grid.grid_elements_pts[i]->max_corner[0],
                grid.grid_elements_pts[i]->min_corner[0],
                grid.grid_elements_pts[i]->min_corner[0],
                grid.grid_elements_pts[i]->max_corner[0],
                grid.grid_elements_pts[i]->min_corner[0],

                grid.grid_elements_pts[i]->min_corner[1],
                grid.grid_elements_pts[i]->min_corner[1],
                grid.grid_elements_pts[i]->min_corner[1],
                grid.grid_elements_pts[i]->max_corner[1],
                grid.grid_elements_pts[i]->max_corner[1],
                grid.grid_elements_pts[i]->max_corner[1],

                grid.grid_elements_pts[i]->min_corner[2],
                grid.grid_elements_pts[i]->max_corner[2],
                grid.grid_elements_pts[i]->max_corner[2],
                grid.grid_elements_pts[i]->min_corner[2],
                grid.grid_elements_pts[i]->max_corner[2],
                grid.grid_elements_pts[i]->max_corner[2],
            };

        vector<double> p1_vals(6, grid.grid_elements_pts[i]->calculated_target_size);

        list_data.insert(list_data.end(), p1_coords.begin(), p1_coords.end());
        list_data.insert(list_data.end(), p1_vals.begin(), p1_vals.end());

        vector<double> p2_coords = {
                grid.grid_elements_pts[i]->min_corner[0],
                grid.grid_elements_pts[i]->max_corner[0],
                grid.grid_elements_pts[i]->max_corner[0],
                grid.grid_elements_pts[i]->min_corner[0],
                grid.grid_elements_pts[i]->max_corner[0],
                grid.grid_elements_pts[i]->max_corner[0],

                grid.grid_elements_pts[i]->min_corner[1],
                grid.grid_elements_pts[i]->min_corner[1],
                grid.grid_elements_pts[i]->min_corner[1],
                grid.grid_elements_pts[i]->max_corner[1],
                grid.grid_elements_pts[i]->max_corner[1],
                grid.grid_elements_pts[i]->max_corner[1],

                grid.grid_elements_pts[i]->min_corner[2],
                grid.grid_elements_pts[i]->max_corner[2],
                grid.grid_elements_pts[i]->min_corner[2],
                grid.grid_elements_pts[i]->min_corner[2],
                grid.grid_elements_pts[i]->max_corner[2],
                grid.grid_elements_pts[i]->min_corner[2],
        };

        vector<double> p2_vals(6, grid.grid_elements_pts[i]->calculated_target_size);

        list_data.insert(list_data.end(), p2_coords.begin(), p2_coords.end());
        list_data.insert(list_data.end(), p2_vals.begin(), p2_vals.end());
    }

    gmsh::view::addListData(view_tag, "SI", grid.grid_elements_pts.size()*2, list_data); // factor of 2 from each grid element cube being split into 2 prisms

    int field_tag = gmsh::model::mesh::field::add("PostView");
    gmsh::model::mesh::field::setNumber(field_tag, "ViewTag", view_tag);

    gmsh::view::write(view_tag, "example_refine_" + to_string(grid.side_length(0)) + ".pos");

    return {(double) field_tag};
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
//    int uniform_field_tag = add_uniform_target_size_field(target_size);
    int uniform_field_tag = add_step_uniform_target_size_field();
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

void analyse_stability_iteration(int iteration, vector<Element*> element_pts, double target_size, double grid_size_lower, double grid_size_upper)
{
    double total_volume = 0;
    double total_target_size = 0;
    double elements_in_sample_volume = 0;

    double min_size = element_pts[0]->linear_target_size;
    double max_size = element_pts[0]->linear_target_size;

    cout << "Num elements: " << element_pts.size() << endl;

    // process mesh
    for (unsigned i=0; i<element_pts.size(); i++)
    {
        // Since we are using a weird shaped object to mesh, not going to select a sample volume - use all elements in analysis:

        total_target_size += element_pts[i]->linear_target_size;
        total_volume += element_pts[i]->volume;
        elements_in_sample_volume += 1;

        if (element_pts[i]->linear_target_size < min_size) min_size = element_pts[i]->linear_target_size;
        if (element_pts[i]->linear_target_size > max_size) max_size = element_pts[i]->linear_target_size;


//        if (element_pts[i]->midpoint_coords[0] > .8 && element_pts[i]->midpoint_coords[0] < 2.2 &&
//            element_pts[i]->midpoint_coords[1] > .8 && element_pts[i]->midpoint_coords[1] < 2.2 &&
//            element_pts[i]->midpoint_coords[2] > .8 && element_pts[i]->midpoint_coords[2] < 2.2)
//        {
//            total_target_size += element_pts[i]->linear_target_size;
//            total_volume += element_pts[i]->volume;
//            elements_in_sample_volume += 1;
//
//            if (element_pts[i]->linear_target_size < min_size) min_size = element_pts[i]->linear_target_size;
//            if (element_pts[i]->linear_target_size > max_size) max_size = element_pts[i]->linear_target_size;
//            cout << "Here" << endl;
//        }
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

    string filename = "grid_stability_t3_bar_bend_twist_its=";
    filename.append(to_string(target_size));
    filename.append("_c=");
    filename.append(to_string(constant));
    filename.append("_N=");
    filename.append(to_string(N));
    filename.append(".txt");

    fstream f;
    f.open(filename, ios_base::app);
    f << iteration << ',' << av_size << ',' << min_size << ',' << max_size << ',' << grid_size_lower << ',' << grid_size_upper << ',' << av_vol << ',' << total_volume << ',' << elements_in_sample_volume << '\n';
    f.close();
}

void run(double target_size)
{
//    double target_size = 0.1;
    string base_out_filename = "out_test";

    cout << "grid based refinement using bg msh, target_size: " << target_size <<  " C = " << constant << ", N = " << N << endl;

    // Make initial uniform mesh
    gmsh::initialize();
    gmsh::clear();

    make_cube(3);

    int uniform_field_tag = add_uniform_target_size_field(target_size);
    mesh_with_bg_fields({(double)uniform_field_tag});

    // Extract mesh data
    vector<Element> elements = get_elements();

    for (unsigned i=0; i<elements.size(); i++)
    {
        elements[i].calculate_midpoint();
        elements[i].calculate_volume();
        elements[i].calculate_target_size_from_volume();
    }

    vector<Element*> elements_pts;
    for (unsigned i=0; i<elements.size(); i++)
    {
        elements_pts.push_back(&elements[i]);
    }

    gmsh::write("t.msh");

    // Once mesh data has been extracted, GMSH is finished with for this iteration.
    gmsh::clear();
    gmsh::finalize();

    cout << "Number of elements: " << elements.size() << endl;

    Grid grid(elements_pts);
    vector<double> initial_min_corner = grid.min_corner;
    vector<double> initial_max_corner = grid.max_corner;

    unsigned n = 10;
    double dx = (initial_max_corner[0] - initial_min_corner[0]) / double(n);


    // shrink grid cube to only include a subset of elements
    for (unsigned iter=0; iter<n; iter++)
    {
        // extract relevant elements
        vector<Element*> new_grid_elements;
        vector<double> new_min_corner(3);
        vector<double> new_max_corner(3);

        double factor = dx * iter;

        for (unsigned i=0; i<3; i++)
        {
            new_min_corner[i] = initial_min_corner[i] + factor;
            new_max_corner[i] = initial_max_corner[i] - factor;
        }

        for (unsigned i=0; i<elements_pts.size(); i++)
        {
            if ((elements_pts[i]->midpoint_coords[0] > new_min_corner[0] && elements_pts[i]->midpoint_coords[0] < new_max_corner[0])
                && (elements_pts[i]->midpoint_coords[1] > new_min_corner[1] && elements_pts[i]->midpoint_coords[1] < new_max_corner[1])
                && (elements_pts[i]->midpoint_coords[2] > new_min_corner[2] && elements_pts[i]->midpoint_coords[2] < new_max_corner[2]))
            {
                new_grid_elements.push_back(elements_pts[i]);
            }
        }

        Grid grid(new_grid_elements);

        gmsh::initialize();
        gmsh::clear();
        add_bg_fields_for_grid(grid);
        gmsh::clear();
        gmsh::finalize();
    }
}

void append_to_file(double target_size, double constant_value)
{
    string filename = "out.txt";
    fstream f;
    f.open(filename, ios_base::app);
    f << target_size << "," << constant_value << "\n";
    f.close();
}

int main()
{
    double target_size = .1;
    run(target_size);

    return 0;
}


//  TODO: Add saving mesh every 10 or 5 iterations so that it can be resumed and analysed at a later date.
