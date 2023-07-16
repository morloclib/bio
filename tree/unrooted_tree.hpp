// requirements:
//   * 0 is the root node
//   * given N == number of nodes and L == number of leafs 
//       indices for nodes range from 0 to (N-1)
//       indices for leafs range from N to N+L-1

#ifndef __MORLOC_BIO_TREE2_HPP__
#define __MORLOC_BIO_TREE2_HPP__

#include <utility>
#include <vector>
#include <variant>
#include <iostream>

#include <mlc.hpp>
#include "tree_types.hpp"


template <typename Node, typename Edge, typename Leaf>
int mlc_unrooted_count_nodes(const UnrootedTree<Node, Edge, Leaf>& tree) {
    int count = 0;
    for (const auto& vertex : tree.verts) {
        if (std::holds_alternative<Node>(vertex)) {
            ++count;
        }
    }
    return count;
}

template <typename Node, typename Edge, typename Leaf>
int mlc_unrooted_count_leafs(const UnrootedTree<Node, Edge, Leaf>& tree) {
    int count = 0;
    for (const auto& vertex : tree.verts) {
        if (std::holds_alternative<Leaf>(vertex)) {
            ++count;
        }
    }
    return count;
}

template <typename Node, typename Edge, typename Leaf>
int mlc_unrooted_count_edges(const UnrootedTree<Node, Edge, Leaf>& tree) {
    int count = 0;
    for (const auto& edges : tree.out) {
        count += edges.size();
    }
    return count;
}

template <typename Node, typename Edge, typename Leaf>
UnrootedTree<Node, Edge, Leaf> mlc_pack_utree(std::tuple<std::vector<Node>, std::vector<std::tuple<int, int, Edge>>, std::vector<Leaf>> pack) {
    const auto& [nodes, edges, leafs] = pack;

    UnrootedTree<Node, Edge, Leaf> tree;

    // Create verts vector from nodes and leafs
    for (const auto& node : nodes) {
        tree.verts.push_back(node);
    }
    for (const auto& leaf : leafs) {
        tree.verts.push_back(leaf);
    }

    // Initialize out
    tree.out.resize(tree.verts.size());

    // Populate out
    for (const auto& edge : edges) {
        int a_index = std::get<0>(edge);
        int b_index = std::get<1>(edge);
        const Edge& edge_value = std::get<2>(edge);

        tree.out[a_index].push_back(std::make_tuple(b_index, edge_value));
        tree.out[b_index].push_back(std::make_tuple(a_index, edge_value));
    }

    return tree;
}

template <typename Node, typename Edge, typename Leaf>
std::tuple<std::vector<Node>, std::vector<std::tuple<int, int, Edge>>, std::vector<Leaf>>
mlc_unpack_utree(const UnrootedTree<Node, Edge, Leaf> &tree)
{
    std::vector<Node> nodes;
    std::vector<Leaf> leafs;
    std::vector<std::tuple<int, int, Edge>> edges;

    for (const auto& vertex : tree.verts) {
        if (std::holds_alternative<Node>(vertex)) {
            nodes.push_back(std::get<Node>(vertex));
        } else {
            leafs.push_back(std::get<Leaf>(vertex));
        }
    }

    for (int i = 0; i < tree.out.size(); ++i) {
        for (const auto& b : tree.out[i]) {
            edges.push_back(std::make_tuple(i, b, std::get<1>(tree.out[i])));
        }
    }

    return {nodes, edges, leafs};
}


// Applies a function to each leaf in the tree and returns the new tree.
template <typename Node, typename Edge, typename Leaf, typename NewLeaf>
UnrootedTree<Node, Edge, NewLeaf> mlc_unrooted_mapLeaf(
    std::function<NewLeaf(Leaf)> func, 
    const UnrootedTree<Node, Edge, Leaf>& tree
) {
    UnrootedTree<Node, Edge, NewLeaf> newTree;
    newTree.out = tree.out;
    newTree.verts.resize(tree.verts.size());

    for(int i=0; i<tree.verts.size(); i++){
        if(std::holds_alternative<Leaf>(tree.verts[i])){
            Leaf oldLeaf = std::get<Leaf>(tree.verts[i]);
            NewLeaf newLeaf = func(oldLeaf);
            newTree.verts[i] = newLeaf;
        }
        else if(std::holds_alternative<Node>(tree.verts[i])){
            Node oldNode = std::get<Node>(tree.verts[i]);
            newTree.verts[i] = oldNode;
        }
    }

    return newTree;
}

// Applies a function to each leaf in the tree and returns the new tree.
template <typename Node, typename Edge, typename Leaf, typename NewNode>
UnrootedTree<NewNode, Edge, Leaf> mlc_unrooted_mapNode(
    std::function<NewNode(Node)> func, 
    const UnrootedTree<Node, Edge, Leaf>& tree
) {
    UnrootedTree<NewNode, Edge, Leaf> newTree;
    newTree.out = tree.out;
    newTree.verts.resize(tree.verts.size());

    for(int i=0; i<tree.verts.size(); i++){
        if(std::holds_alternative<Node>(tree.verts[i])){
            Node oldNode = std::get<Node>(tree.verts[i]);
            NewNode newNode = func(oldNode);
            newTree.verts[i] = newNode;
        }
        else if(std::holds_alternative<Leaf>(tree.verts[i])){
            Leaf oldLeaf = std::get<Leaf>(tree.verts[i]);
            newTree.verts[i] = oldLeaf;
        }
    }

    return newTree;
}

// Applies a function to each edge in the tree and returns the new tree.
template <typename Node, typename Edge, typename Leaf, typename NewEdge>
UnrootedTree<Node, NewEdge, Leaf> mlc_unrooted_mapEdge(
    std::function<NewEdge(Edge)> func, 
    const UnrootedTree<Node, Edge, Leaf>& tree
) {
    UnrootedTree<NewNode, Edge, Leaf> newTree;
    newTree.vert = tree.vert;

    for(int i=0; i<tree.out.size(); i++){
        // get child index
        int childIndex = std::get<0>(tree.out[i]);

        // get child original edge value
        Edge newEdge = func(std::get<1>(tree.out[i])); 

        newTree.out.push_back(std::make_tuple(childIndex, newEdge));
    }

    return newTree;
}

template<typename Node, typename Edge, typename Leaf, typename B>
UnrootedTree<Node, Edge, Leaf> mlc_unrooted_treeBy(
    std::function<UnrootedTree<Node, Edge, int>(std::vector<B>)> buildTree,
    std::vector<std::tuple<Leaf, B>> xs
) {
    std::vector<B> buildTreeArgs;
    for (auto& pair : xs) {
        buildTreeArgs.push_back(std::get<1>(pair));
    }

    UnrootedTree<Node, Edge, int> indexTree = buildTree(buildTreeArgs);

    UnrootedTree<Node, Edge, Leaf> finalTree = mlc_unrooted_mapLeaf(
        [&](int index) { return std::get<0>(xs[index]); },
        indexTree
    );

    return finalTree;
}


template <typename Node, typename Edge, typename Leaf>
RootedTree<Node, Edge, Leaf> dfsRooting(
    UnrootedTree<Node, Edge, Leaf>& unrootedTree,
    int currentNodeIndex,
    int parentIndex
) {
    RootedTree<Node, Edge, Leaf> rootedTree;
    rootedTree.data = std::get<Node>(unrootedTree.verts[currentNodeIndex]);

    // Visit all neighbors of the current node
    for (int neighborIndex : unrootedTree.out[currentNodeIndex]) {
        // Skip the parent node to prevent cycles
        if (neighborIndex == parentIndex) {
            continue;
        }

        // If the neighbor is a node, recursively root it
        if (std::holds_alternative<Node>(unrootedTree.verts[neighborIndex])) {
            RootedTree<Node, Edge, Leaf> subtree = dfsRooting(unrootedTree, neighborIndex, currentNodeIndex);
            rootedTree.children.push_back(subtree);
            rootedTree.edges.push_back(std::get<1>(unrootedTree.edges[neighborIndex]));
        }

        // If the neighbor is a leaf, add it directly to the children
        else if (std::holds_alternative<Leaf>(unrootedTree.verts[neighborIndex])) {
            Leaf leaf = std::get<Leaf>(unrootedTree.verts[neighborIndex]);
            rootedTree.children.push_back(leaf);
            rootedTree.edges.push_back(unrootedTree.edges[neighborIndex]);
        }
    }

    return rootedTree;
}

template <typename Node, typename Edge, typename Leaf>
RootedTree<Node, Edge, Leaf> mlc_unrooted_root_at(
    UnrootedTree<Node, Edge, Leaf>& unrootedTree, 
    int rootIndex
) {
    return dfsRooting(unrootedTree, rootIndex, -1);
}




// mutating helper function to add an edge
template <typename Node, typename Leaf>
void _add_edge(int a, int b, Edge e, UnrootedTree<Node, Edge, Leaf>& unrootedTree){
    unrootedTree.out[a].push_back(std::make_tuple(b, e));
    unrootedTree.out[b].push_back(std::make_tuple(a, e));
}

// mutating helper function to split an edge and return the new index
template <typename Node, typename Edge, typename Leaf>
int _split_edge(int a, int b, Edge ea, Edge eb, Node n, UnrootedTree<Node, Edge, Leaf>& unrootedTree){

    // this will be the index of the new node
    int newIndex = unrootedTree.out.size();

    // push the new node
    unrootedTree.verts.push_back(n);

    // link from the new node to a
    _add_edge(newIndex, a, ea, unrootedTree);

    // link from the new node to b
    _add_edge(newIndex, b, eb, unrootedTree);

    return newIndex;
}


template <typename Node, typename Leaf>
std::tuple<std::vector<int>, std::vector<double>, double> _find_distal_leaf_r(const UnrootedTree<Node, double, Leaf>& unrootedTree, int nodeIndex, int parentIndex){
    std::vector<int> bestPath;
    std::vector<double> distances;
    double maxDist = 0.0;
    double bestChildEdge = -1.0;

    for(int i = 0; i < unrootedTree.out.size(); i++){
        const auto& [childIndex, childEdge] = unrootedTree.out[i];

        // Should we infinite loop? Let's not
        if (parentIndex == childIndex){
            continue;
        }

        // Shall we plumb the depths? Hell yes
        const auto& [childPath, childDistances, childMaxDist] = _find_distal_leaf_r(unrootedTree, childIndex, nodeIndex);

        // Are we happy now? We'll see
        if ((childMaxDist + childEdge) > maxDist){
            // So no, hell yeah, replace it
            maxDist = childMaxDist + childEdge;
            distances = childDistances;
            bestPath = childPath;
            bestChildEdge = childEdge;
        }
    }

    if (bestChildEdge > 0){
        distances.push_back(bestChildEdge);
    }

    bestPath.push_back(nodeIndex);

    // And pass it up
    return std::make_tuple(bestPath, distances, maxDist);
}

template <typename Node, typename Leaf>
std::tuple<std::vector<int>, std::vector<double>, double> _find_distal_leaf(const UnrootedTree<Node, double, Leaf>& unrootedTree, int nodeIndex){
    return _find_distal_leaf_r(unrootedTree, nodeIndex, -1);
}

// For a tree, finding the longest path can be done by finding the most distant
// leaf from a randome node and then finding the longest path starting from that node.
template <typename Node, typename Leaf>
std::tuple<std::vector<int>, std::vector<double>, double> findLongestPath(const UnrootedTree<Node, double, Leaf>& unrootedTree){
    // Find the longest path by summed edge length from a random start. I use 0
    // here, but it could be anything.
    const auto& [path_a, path_length_a] = _find_distal_leaf(unrootedTree, 0);

    // Find the longest path from this leaf
    return _find_distal_leaf(unrootedTree, path_a);
}

template <typename Node, typename Leaf>
RootedTree<Node, double, Leaf> mlc_unrooted_midpoint(UnrootedTree<Node, double, Leaf>& unrootedTree, Node rootNode) {

    // Find the longest path (by total edge length) in the unrooted tree.
    // The returned path should be a list of vertices in the order they are visited in the path.
    const auto& [longestPath, distances, totalLength] = findLongestPath(unrootedTree);

    float pathTraversed = 0.0;
    int rootIndex;
    int a = longestPath[0];
    int b;

    for(int i = 0; i < distances.size(); i++){
        float edge = distances[i];
        b = longestPath[i+1];
        if (pathTraversed + edge > (totalLength / 2)){
            // Break the current branch at the midpoint
            float ea = (totalLength / 2) - pathTraversed;
            float eb = (pathTraversed + edge) - (totalLength / 2);

            // Set the new root index (we are just adding one new node
            int rootIndex = _split_edge(a, b, ea, eb, rootNode, unrootedTree);

            // Add a root at the midpoint of the longest path in the unrooted tree.
            // This will return a rooted tree.
            return mlc_unrooted_root_at(unrootedTree, rootIndex);
        }
        a = b;
    }
}


#endif
