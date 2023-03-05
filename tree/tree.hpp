// requirements:
//   * 0 is the root node
//   * given N == number of nodes and L == number of leafs 
//       indices for nodes range from 0 to (N-1)
//       indices for leafs range from N to N+L-1

#ifndef __TREE_HPP__
#define __TREE_HPP__

#include <utility>
#include <vector>
#include <variant>

template <typename Node, typename Edge, typename Leaf>
struct Tree {
    Node data;
    std::vector<std::variant<Tree, Leaf>> children;
    std::vector<Edge> edges;
};

template <typename Node, typename Edge, typename Leaf>
int mlc_count_nodes(const Tree<Node, Edge, Leaf>& tree) {
    int count = 1;
    for (const auto& child : tree.children) {
        if (std::holds_alternative<Tree<Node, Edge, Leaf>>(child)) {
            count += mlc_count_nodes(std::get<Tree<Node, Edge, Leaf>>(child));
        }
    }
    return count;
}

template <typename Node, typename Edge, typename Leaf>
Tree<Node, Edge, Leaf> mlc_pack_tree_r(
        size_t node_index,
        const std::vector<Node>& nodes,
        const std::vector<Leaf>& leafs,
        const std::vector<std::vector<Edge>>& edges,
        const std::vector<std::vector<int>>& child_indices) {
    std::vector<std::variant<Tree<Node, Edge, Leaf>, Leaf>> children;

    for (int child_index : child_indices[node_index]) {
        if (child_index < nodes.size()) {
            // If it is a node
            Tree<Node, Edge, Leaf> child = mlc_pack_tree_r(
                child_index, nodes, leafs, edges, child_indices);
            children.push_back(child);
        } else {
            // If it is a leaf
            Leaf leaf = leafs[child_index - nodes.size()];
            children.push_back(leaf);
        }
    }

    return Tree<Node, Edge, Leaf>{nodes[node_index], children, edges[node_index]};
}

template <typename Node, typename Edge, typename Leaf>
Tree<Node, Edge, Leaf> mlc_pack_tree(const std::tuple<std::vector<Node>,
                           std::vector<std::tuple<int, int, Edge>>,
                           std::vector<Leaf>>& pack) {
    const auto& [nodes, edges, leafs] = pack;

    // Store the children and edges for every node (these include leaf and node children)
    std::vector<std::vector<int>> child_indices(nodes.size());
    std::vector<std::vector<Edge>> node_edges(nodes.size());
    for (const auto& edge : edges) {
        int parent_index = std::get<0>(edge);
        int child_index = std::get<1>(edge);
        const Edge& edge_value = std::get<2>(edge);
        node_edges[parent_index].push_back(edge_value);
        child_indices[parent_index].push_back(child_index);
    }

    // The 0 index is root
    return mlc_pack_tree_r(0, nodes, leafs, node_edges, child_indices);
}


template <typename Node, typename Edge, typename Leaf>
std::tuple<std::vector<Node>, std::vector<std::tuple<int, int, Edge>>, std::vector<Leaf>>
mlc_unpack_tree(const Tree<Node, Edge, Leaf>& tree)
{
    std::vector<Node> nodes;
    std::vector<Leaf> leafs;
    std::vector<std::tuple<int, int, Edge>> edges;

    int number_of_nodes = mlc_count_nodes(tree);
    mlc_unpack_tree_r(tree, nodes, leafs, edges, number_of_nodes, 0);
    return {nodes, edges, leafs};
}


template <typename Node, typename Edge, typename Leaf>
void mlc_unpack_tree_r(const Tree<Node, Edge, Leaf>& tree, std::vector<Node>& nodes,
                       std::vector<Leaf>& leafs, std::vector<std::tuple<int, int, Edge>>& edges, int number_of_nodes, int index)
{
    nodes.push_back(tree.data);
    for (const auto& child : tree.children) {
        if (std::holds_alternative<Tree<Node, Edge, Leaf>>(child)) {
            const auto& child_tree = std::get<Tree<Node, Edge, Leaf>>(child);
            edges.push_back(std::make_tuple(index, nodes.size(), child_tree.edges.front()));
            mlc_unpack_tree_r(child_tree, nodes, leafs, edges, number_of_nodes, nodes.size());
        } else {
            leafs.push_back(std::get<Leaf>(child));
            edges.push_back(std::make_tuple(index, leafs.size() + number_of_nodes - 1, tree.edges.front()));
        }
    }
}

#endif
