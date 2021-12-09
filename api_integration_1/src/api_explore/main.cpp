//
// Created by Ben Gavan on 21/07/2021.
//


#include <iostream>
#include <string>
#include "../../libs/gmsh.h"
//#include <gmsh.h>

using namespace std;

class Node {
public:
    size_t tag;
    std::vector<double> coord;
    std::vector<double> parametricCoord;

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

// get all elementary entities in the model
std::vector <std::pair<int, int>> getModelEntities()
{
    std::vector <std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);
    return entities;
}

void getMeshElements()
{
    using namespace std;
    vector<int> elemTypes;
    vector <vector<size_t>> elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags);
    for (auto e : elemTypes) {
        cout << e << '\n';
    }
    cout << "elemNodeTags" << endl;
    for (vector<size_t> v : elemNodeTags) {
        cout << "size: " << v.size() << "\n";
        for (size_t s : v) {
            cout << s << '\n';
        }
    }
}

std::vector<std::size_t> getTetrahedronMeshElementTags()
{
    std::vector<std::size_t> elementTags;
    std::vector<std::size_t> nodeTags;
    // 4 = 4-node tetrahedron
    gmsh::model::mesh::getElementsByType(4, elementTags, nodeTags);
    return elementTags;
}

Node getMeshNode(const size_t nodeTag)
{
    std::vector<double> coord;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNode(nodeTag, coord, parametricCoord);
    return Node(nodeTag, coord, parametricCoord);
}

vector<Node> getMeshNodes()
{
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords;
    std::vector<double> parametricCoords;

    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoords);

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


void mine() {
    gmsh::initialize();

    string filename = "../../out/cube_withcallback_1.msh";

    gmsh::open(filename);

    cout << "opened" << '\n';

//    vector <pair<int, int>> entities = getModelEntities();
//
//    cout << entities.size() << endl;

//    cout << "Get mesh elements" << endl;
//    getMeshElements();

    // Get Tetrahedron mesh element tags and the tags of the constituent nodes
    std::vector<std::size_t> elementTags;
    std::vector<std::size_t> nodeTags;
    // 4 = 4-node tetrahedron
    gmsh::model::mesh::getElementsByType(4, elementTags, nodeTags);

    int num_elements = elementTags.size();


    // can use 'getNode(...)' to get specific node coords & parametric coords
    // but "This function relies on an internal cache (a vector in
    // case of dense node numbering, a map otherwise); for large meshes accessing
    // nodes in bulk is often preferable."

    // access all nodes again (doing in bulk (due to internal cache))
    cout << "getMeshNodes()\n";

    vector<Node> nodes = getMeshNodes();

    cout << "here" << endl;

    cout << "nodes.size(): " << nodes.size() << endl;

    cout << "Nodes tags: ";
    for (Node n : nodes) {
        cout << n.tag << ": ";
        n.printCoords();
        cout << '\n';
    }
    cout << endl;



//    vector<vector<double>>

    gmsh::finalize();
}

void old() {
    gmsh::initialize();

    string filename = "../../out/cube_withcallback_1.msh";

    gmsh::open(filename);

    // get all elementary entities in the model
    std::vector <std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);

    for (unsigned int i = 0; i < entities.size(); i++) {
        // get the mesh nodes for each elementary entity
        std::vector <std::size_t> nodeTags;
        std::vector<double> nodeCoords, nodeParams;
        int dim = entities[i].first, tag = entities[i].second;
        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);

        // get the mesh elements for each elementary entity
        std::vector<int> elemTypes;
        std::vector <std::vector<std::size_t>> elemTags, elemNodeTags;
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

        // report some statistics
        int numElem = 0;
        for (unsigned int j = 0; j < elemTags.size(); j++)
            numElem += elemTags[j].size();
        std::string type;
        gmsh::model::getType(dim, tag, type);
        std::cout << nodeTags.size() << " mesh nodes and " << numElem
                  << " mesh elements on entity (" << dim << "," << tag
                  << ") of type " << type << "\n";
        std::vector<int> partitions;
        gmsh::model::getPartitions(dim, tag, partitions);
        if (partitions.size()) {
            std::cout << " - Partition tag(s):";
            for (unsigned int j = 0; j < partitions.size(); j++)
                std::cout << " " << partitions[j];
            int parentDim, parentTag;
            gmsh::model::getParent(dim, tag, parentDim, parentTag);
            std::cout << " - parent entity (" << parentDim << "," << parentTag
                      << ")\n";
        }
        for (unsigned int j = 0; j < elemTypes.size(); j++) {
            std::string name;
            int d, order, numv, numpv;
            std::vector<double> param;
            gmsh::model::mesh::getElementProperties(elemTypes[j], name, d, order,
                                                    numv, param, numpv);
            std::cout << " - Element type: " << name << ", order " << order << "\n";
            std::cout << "   with " << numv << " nodes in param coord: (";
            for (unsigned int k = 0; k < param.size(); k++)
                std::cout << param[k] << " ";
            std::cout << ")\n";
        }
    }

    gmsh::finalize();
}

int main(int argc, char **argv) {
    mine();
    return 0;
}

