#ifndef __MORLOC_BIO_DIRECTED_TREE_HPP__
#define __MORLOC_BIO_DIRECTED_TREE_HPP__

#include <utility>
#include <vector>
#include <variant>
#include <iostream>

#include <mlc.hpp>

template <typename Node, typename Edge, typename Leaf>
struct DirectedTree {
    std::vector<std::vector<std::tuple<int, Edge>>> edge;
    std::vector<std::variant<Node,Leaf>> verts;
};

template <typename Node, typename Edge, typename Leaf>
int mlc_directed_count_nodes(const DirectedTree<Node, Edge, Leaf>& tree) {
    int count = 0;
    for (const auto& vertex : tree.verts) {
        if (std::holds_alternative<Node>(vertex)) {
            ++count;
        }
    }
    return count;
}

template <typename Node, typename Edge, typename Leaf>
int mlc_directed_count_leafs(const DirectedTree<Node, Edge, Leaf>& tree) {
    int count = 0;
    for (const auto& vertex : tree.verts) {
        if (std::holds_alternative<Leaf>(vertex)) {
            ++count;
        }
    }
    return count;
}

template <typename Node, typename Edge, typename Leaf>
DirectedTree<Node, Edge, Leaf> mlc_pack_directed_tree(std::tuple<std::vector<Node>, std::vector<std::tuple<int, int, Edge>>, std::vector<Leaf>> pack) {
    const auto& [nodes, edges, leafs] = pack;

    DirectedTree<Node, Edge, Leaf> tree;

    // Create verts vector from nodes and leafs
    for (const auto& node : nodes) {
        tree.verts.push_back(node);
    }
    for (const auto& leaf : leafs) {
        tree.verts.push_back(leaf);
    }

    // Initialize edges vector with empty vectors
    tree.edges.resize(tree.verts.size());

    // Populate edges vector using input edges
    for (const auto& edge : edges) {
        int parent_index = std::get<0>(edge);
        int child_index = std::get<1>(edge);
        const Edge& edge_value = std::get<2>(edge);
        tree.edges[parent_index].push_back(std::make_tuple(child_index, edge_value));
    }

    return tree;
}

template <typename Node, typename Edge, typename Leaf>
std::tuple<std::vector<Node>, std::vector<std::tuple<int, Edge>>, std::vector<Leaf>>
mlc_unpack_directed_tree(const DirectedTree<Node, Edge, Leaf> &tree)
{
    std::vector<Node> nodes;
    std::vector<Leaf> leafs;
    std::vector<std::tuple<int, Edge>> edges;

    for (const auto& vertex : tree.verts) {
        if (std::holds_alternative<Node>(vertex)) {
            nodes.push_back(std::get<Node>(vertex));
        } else {
            leafs.push_back(std::get<Leaf>(vertex));
        }
    }

    for (const auto& edge_list : tree.edges) {
        for (const auto& edge : edge_list) {
            edges.push_back(edge);
        }
    }

    return {nodes, edges, leafs};
}

// Applies a function to each leaf in the tree and returns the new tree.
// mapLeaf :: (l -> l') -> Tree n e l -> Tree n e l'
template <typename Node, typename Edge, typename Leaf, typename NewLeaf>
DirectedTree<Node, Edge, NewLeaf> mlc_directed_mapLeaf(
    std::function<NewLeaf(Leaf)> func, 
    DirectedTree<Node, Edge, Leaf>& tree
) {
    DirectedTree<Node, Edge, NewLeaf> newTree;
    newTree.edges = tree.edges;
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
// mapNode :: (l -> l') -> Tree n e l -> Tree n e l'
template <typename Node, typename Edge, typename Leaf, typename NewNode>
DirectedTree<Node, Edge, NewNode> mlc_directed_mapNode(
    std::function<NewNode(Node)> func, 
    DirectedTree<Node, Edge, Leaf>& tree
) {
    DirectedTree<Node, Edge, NewLeaf> newTree;
    newTree.edges = tree.edges;
    newTree.verts.resize(tree.verts.size());

    for(int i=0; i<tree.verts.size(); i++){
        if(std::holds_alternative<Node>(tree.verts[i])){
            Leaf oldNode = std::get<Node>(tree.verts[i]);
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

#ifndef
