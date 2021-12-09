//
// Created by Ben Gavan on 17/07/2021.
//

#include <vector>
#include <iostream>
#include <string>
#include <tuple>
#include <cmath>

using namespace std;

class Node {
public:
    double x, y, z;
    double val;
    Node(const double &x, const double &y, const double &z, const double &val);
    void print();
};

Node::Node(const double &x, const double &y, const double &z, const double &val) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->val = val;
    return;
}

void Node::print() {
    printf("(%f, %f, %f), v = %f\n", x, y, z, val);
}

class Mesh {
private:
public:
    vector<Node> nodes;
    Mesh(vector<Node> nodes);
    Mesh newRandomMesh();
    void printNodes();
};

Mesh::Mesh(vector<Node> nodes) {
    this->nodes = nodes;
}

void Mesh::printNodes() {
    for (auto n : this->nodes) {
        n.print();
    }
}

Mesh newRandomMesh() {
    Node n1 (0,0,0,1);
    Node n2 (1,2,3,2);
    Node n3 (-1,-2,-3,1);
    vector<Node> nodes = {n1, n2, n3};
    Mesh m (nodes);
    return m;
}

class Grid {
private:
    Mesh mesh = Mesh({});
    vector<double> min;
    vector<double> max;
    tuple<vector<double>, vector<double>> minMax();
public:
    Grid(Mesh m);
    void printMesh();
};

// Finds the min and max coordinates to define the total size of the cuboid which is then separated into cells.
tuple<vector<double>, vector<double>> Grid::minMax() {
    Node* n = &mesh.nodes[0];

    vector<double> min = {n->x, n->y, n->z};
    vector<double> max = {n->x, n->y, n->z};

    for (unsigned i = 1; i < mesh.nodes.size(); i++) {
        n = &mesh.nodes[i];
        if (n->x < min[0]) min[0] = n->x;
        if (n->x > max[0]) max[0] = n->x;

        if (n->y < min[1]) min[1] = n->y;
        if (n->y > max[1]) max[1] = n->y;

        if (n->z < min[2]) min[2] = n->z;
        if (n->z > max[2]) max[2] = n->z;
    }
    return make_tuple(min, max);
}

// Makes the grid for node finding
Grid::Grid(Mesh m) {
    mesh = m;

    tie(min, max) = this->minMax();

    // Calculate the overall size of the total grid
    float r = std::pow(m.nodes.size(), 1/3.); //

    int n = ceil(r); //
    float cf = ceil(r);
    cout << "c = " << c << endl;
    cout << "cf = " << cf << endl;

    float l_x = (max[0] - min[0]) / r; // length of x side of each cell
    float l_y = (max[1] - min[1]) / r;
    float l_z = (max[2] - min[2]) / r;

    cout << r << ", " << l_x << ", " << l_y << ", " << l_z << '\n';

    vector<vector<vector<vector<Node*>>>> grid;

    for (Node n : mesh.nodes) {

    }


}

void Grid::printMesh() {
    mesh.printNodes();
}





int main()
{
    Mesh m = newRandomMesh();
    m.printNodes();

    // split mesh into 3d vector to node finding based on position
    Grid g (m);
//    g.printMesh();
    // get relevant nodes for search point/coord

    // interpolate value
    return 0;
}