// requirements:
//   * 0 is the root node
//   * given N == number of nodes and L == number of leafs 
//       indices for nodes range from 0 to (N-1)
//       indices for leafs range from N to N+L-1

#ifndef __MORLOC_BIO_ROOTED_TREE_HPP__
#define __MORLOC_BIO_ROOTED_TREE_HPP__

#include <utility>
#include <vector>
#include <variant>
#include <iostream>

#include "mlccpptypes/prelude.hpp"
#include "mlccpptypes/bio.hpp"

template <typename Node, typename Edge, typename Leaf>
int mlc_rooted_count_nodes(RootedTree<Node, Edge, Leaf> tree) {
    int count = 1;
    for (const auto& child : tree.children) {
        if (std::holds_alternative<RootedTree<Node, Edge, Leaf>>(child)) {
            count += mlc_rooted_count_nodes(std::get<RootedTree<Node, Edge, Leaf>>(child));
        }
    }
    return count;
}

template <typename Node, typename Edge, typename Leaf>
RootedTree<Node, Edge, Leaf> mlc_rooted_pack_tree_r(
        size_t node_index,
        std::vector<Node> nodes,
        std::vector<Leaf> leafs,
        std::vector<std::vector<Edge>> edges,
        std::vector<std::vector<int>> child_indices) {
    std::vector<std::variant<RootedTree<Node, Edge, Leaf>, Leaf>> children;

    for (int child_index : child_indices[node_index]) {
        if (child_index < nodes.size()) {
            // If it is a node
            RootedTree<Node, Edge, Leaf> child = mlc_rooted_pack_tree_r(
                child_index, nodes, leafs, edges, child_indices);
            children.push_back(child);
        } else {
            // If it is a leaf
            Leaf leaf = leafs[child_index - nodes.size()];
            children.push_back(leaf);
        }
    }

    return RootedTree<Node, Edge, Leaf>{nodes[node_index], children, edges[node_index]};
}

template <typename Node, typename Edge, typename Leaf>
RootedTree<Node, Edge, Leaf> mlc_rooted_pack_tree(std::tuple<std::vector<Node>,
                           std::vector<std::tuple<int, int, Edge>>,
                           std::vector<Leaf>> pack) {
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
    return mlc_rooted_pack_tree_r(0, nodes, leafs, node_edges, child_indices);
}


template <typename Node, typename Edge, typename Leaf>
std::tuple<std::vector<Node>, std::vector<std::tuple<int, int, Edge>>, std::vector<Leaf>>
mlc_rooted_unpack_tree(const RootedTree<Node, Edge, Leaf> &tree)
{
    std::vector<Node> nodes;
    std::vector<Leaf> leafs;
    std::vector<std::tuple<int, int, Edge>> edges;

    int number_of_nodes = mlc_rooted_count_nodes(tree);
    mlc_rooted_unpack_tree_r(tree, nodes, leafs, edges, number_of_nodes, 0);
    return {nodes, edges, leafs};
}


template <typename Node, typename Edge, typename Leaf>
void mlc_rooted_unpack_tree_r(const RootedTree<Node, Edge, Leaf>& tree, std::vector<Node>& nodes,
                       std::vector<Leaf>& leafs, std::vector<std::tuple<int, int, Edge>>& edges, int number_of_nodes, int index)
{
    nodes.push_back(tree.data);
    for (const auto& child : tree.children) {
        if (std::holds_alternative<RootedTree<Node, Edge, Leaf>>(child)) {
            const auto& child_tree = std::get<RootedTree<Node, Edge, Leaf>>(child);
            edges.push_back(std::make_tuple(index, nodes.size(), child_tree.edges.front()));
            mlc_rooted_unpack_tree_r(child_tree, nodes, leafs, edges, number_of_nodes, nodes.size());
        } else {
            leafs.push_back(std::get<Leaf>(child));
            edges.push_back(std::make_tuple(index, leafs.size() + number_of_nodes - 1, tree.edges.front()));
        }
    }
}

// node :: RootedTree n e l -> n
template <typename Node, typename Edge, typename Leaf>
Node mlc_rooted_node(RootedTree<Node, Edge, Leaf> tree){
    return tree.data;
}

// subtrees :: RootedTree n e l -> [(e, RootedTree n e l)]
template <typename Node, typename Edge, typename Leaf>
std::vector<std::tuple<Edge, RootedTree<Node, Edge, Leaf>>> mlc_rooted_subtrees(RootedTree<Node, Edge, Leaf> tree){
    std::vector<std::tuple<Edge, RootedTree<Node, Edge, Leaf>>> subtrees;
    for(std::size_t i = 0; i < tree.children.size(); i++){
        if (std::holds_alternative<RootedTree<Node,Edge,Leaf>>(tree.children[i])) {
            subtrees.push_back(std::make_tuple(
              tree.edges[i],
              std::get<RootedTree<Node,Edge,Leaf>>(tree.children[i])
            ));
        }
    }
    return subtrees;
}

// childLeafs :: RootedTree n e l -> [(e, l)]
template <typename Node, typename Edge, typename Leaf>
std::vector<std::tuple<Edge, Leaf>> mlc_rooted_childLeafs(RootedTree<Node, Edge, Leaf> tree){
    std::vector<std::tuple<Edge, Leaf>> leafs;
    for(std::size_t i = 0; i < tree.children.size(); i++){
        if (std::holds_alternative<Leaf>(tree.children[i])) {
            leafs.push_back(std::make_tuple(
              tree.edges[i],
              std::get<Leaf>(tree.children[i])
            ));
        }
    }
    return leafs;
}

template <typename Node, typename Edge, typename Leaf, typename NewLeaf>
RootedTree<Node, Edge, NewLeaf> mlc_rooted_mapLeaf(
    std::function<NewLeaf(Leaf)> func, 
    RootedTree<Node, Edge, Leaf> tree
) {
    RootedTree<Node, Edge, NewLeaf> newRootedTree;
    newRootedTree.data = tree.data;
    newRootedTree.edges = tree.edges;

    for (const auto& child : tree.children) {
        if (std::holds_alternative<Leaf>(child)) {
            // It's a leaf. Apply the function to it and add it to newRootedTree.
            Leaf oldLeaf = std::get<Leaf>(child);
            NewLeaf newLeaf = func(oldLeaf);
            newRootedTree.children.push_back(newLeaf);
        }
        else if (std::holds_alternative<RootedTree<Node, Edge, Leaf>>(child)) {
            // It's a tree. Recursively call mapLeafs on it and add it to newRootedTree.
            const RootedTree<Node, Edge, Leaf>& oldSubtree = std::get<RootedTree<Node, Edge, Leaf>>(child);
            RootedTree<Node, Edge, NewLeaf> newSubtree = mlc_rooted_mapLeaf(func, oldSubtree);
            newRootedTree.children.push_back(newSubtree);
        }
    }

    return newRootedTree;
}


// The list xs contains pairs of values. The first element in the pair is the
// value that will ultimately be stored in the tree's leafs. The second element
// is data that is used to generate the tree. The list xs is split into a list
// of B's that is passed to the buildTree function. The output tree contains
// integer indices as leaves. These indices are used to lookup the appropriate
// leaf in xs and replace the list in the tree.
//
// treeBy :: ([b] -> RootedTree n e Int) -> [(l, b)] -> RootedTree n e l
template<typename Node, typename Edge, typename Leaf, typename B>
RootedTree<Node, Edge, Leaf> mlc_rooted_treeBy(
    std::function<RootedTree<Node, Edge, int>(std::vector<B>)> buildTree,
    std::vector<std::tuple<Leaf, B>> xs
) {
    // Extract the second element of each pair in xs to pass to buildTree.
    std::vector<B> buildTreeArgs;
    for (auto& pair : xs) {
        buildTreeArgs.push_back(std::get<1>(pair));
    }

    // Build the initial tree.
    RootedTree<Node, Edge, int> indexTree = buildTree(buildTreeArgs);

    // Replace indices with their corresponding leafs
    RootedTree<Node, Edge, Leaf> finalTree = mlc_rooted_mapLeaf(
        static_cast<std::function<Leaf(int)>>(
            [&](int index) {
               return std::get<0>(xs[index]);
            }
        ),
        indexTree
    );

    return finalTree;
}


// push :: (n -> n')
//      -> (n' -> e -> n -> (e', n'))
//      -> (n' -> e -> l -> (e', l'))
//      -> RootedTree n e l
//      -> RootedTree n' e' l'
template<typename Node, typename Edge, typename Leaf, typename NodePrime, typename EdgePrime, typename LeafPrime>
RootedTree<NodePrime, EdgePrime, LeafPrime> mlc_rooted_push(
    std::function<NodePrime(Node)> handleRoot,
    std::function<std::tuple<EdgePrime, NodePrime>(NodePrime, Edge, Node)> alterChildNode,
    std::function<std::tuple<EdgePrime, LeafPrime>(NodePrime, Edge, Leaf)> alterLeaf,
    RootedTree<Node, Edge, Leaf> oldRoot
) {
    NodePrime newRootNode = handleRoot(oldRoot.data);
    return mlc_rooted_push_r(alterChildNode, alterLeaf, newRootNode, oldRoot);
}

// push_r :: (n' -> e -> n -> (e', n'))
//        -> (n' -> e -> l -> (e', l'))
//        -> n'
//        -> RootedTree n e l
//        -> RootedTree n' e' l'
template<typename Node, typename Edge, typename Leaf, typename NodePrime, typename EdgePrime, typename LeafPrime>
RootedTree<NodePrime, EdgePrime, LeafPrime> mlc_rooted_push_r(
    std::function<std::tuple<EdgePrime, NodePrime>(NodePrime, Edge, Node)> alterChildNode,
    std::function<std::tuple<EdgePrime, LeafPrime>(NodePrime, Edge, Leaf)> alterLeaf,
    NodePrime newNode,
    RootedTree<Node, Edge, Leaf> oldTree
) {
    RootedTree<NodePrime, EdgePrime, LeafPrime> newTree;
    newTree.data = newNode;
    for(std::size_t i = 0; i < oldTree.children.size(); i++){
        auto child = oldTree.children[i];
        auto edge = oldTree.edges[i];
        if (std::holds_alternative<Leaf>(child)) {
            Leaf oldLeaf = std::get<Leaf>(child);
            auto newEdgeAndLeaf = alterLeaf(newNode, edge, oldLeaf);
            newTree.edges.push_back(std::get<0>(newEdgeAndLeaf));
            newTree.children.push_back(std::get<1>(newEdgeAndLeaf));
        }
        else if (std::holds_alternative<RootedTree<Node, Edge, Leaf>>(child)) {
            RootedTree<Node, Edge, Leaf> oldSubtree = std::get<RootedTree<Node, Edge, Leaf>>(child);
            auto newEdgeAndNode = alterChildNode(newNode, edge, oldSubtree.data);
            RootedTree<NodePrime, EdgePrime, LeafPrime> newSubtree = mlc_rooted_push_r(alterChildNode, alterLeaf, std::get<1>(newEdgeAndNode), oldSubtree);
            newTree.edges.push_back(std::get<0>(newEdgeAndNode));
            newTree.children.push_back(newSubtree);
        }
    }
    return newTree;
}


// pull :: (l -> n')
//      -> (n -> e -> n' -> e')
//      -> (n -> [(e', n')] -> n')
//      -> RootedTree n e l
//      -> RootedTree n' e' l
template<typename Leaf, typename NodePrime, typename Node, typename Edge, typename EdgePrime>
RootedTree<NodePrime, EdgePrime, Leaf> mlc_rooted_pull(
    std::function<NodePrime(Leaf)> handleLeaf,
    std::function<EdgePrime(Node, Edge, NodePrime)> updateEdge,
    std::function<NodePrime(Node, std::vector<std::tuple<EdgePrime, NodePrime>>)> updateNode,
    RootedTree<Node, Edge, Leaf> initialTree 
) {
    RootedTree<NodePrime, EdgePrime, Leaf> newTree;

    std::vector<std::tuple<Edge, NodePrime>> links;

    for (std::size_t i = 0; i < initialTree.children.size(); i++){
        auto child = initialTree.children[i];
        auto edge = initialTree.edges[i];
        if (std::holds_alternative<Leaf>(child)) {
            // It's a leaf. Do not change the leaf, keep it in the new tree.
            // But from the leaf, synthesize a node value that will be
            // used to synthesize parent nodes.
            Leaf oldLeaf = std::get<Leaf>(child);
            NodePrime synthesizedNode = handleLeaf(oldLeaf);
            newTree.children.push_back(oldLeaf);
            links.push_back(std::make_tuple(edge, synthesizedNode));
        }
        else if (std::holds_alternative<RootedTree<Node, Edge, Leaf>>(child)) {
            // It's a tree. Recursively call mlc_rooted_pull on it and add it to newTree.
            RootedTree<Node, Edge, Leaf> oldSubtree = std::get<RootedTree<Node, Edge, Leaf>>(child);
            RootedTree<NodePrime, EdgePrime, Leaf> newSubtree = mlc_rooted_pull(handleLeaf, updateEdge, updateNode, oldSubtree);
            newTree.children.push_back(newSubtree);
            links.push_back(std::make_tuple(edge, newSubtree.data));
        }
    }

    std::vector<std::tuple<EdgePrime, NodePrime>> updatedLinks;

    for (auto link : links){
        // synthesize new edge
        EdgePrime newEdge = updateEdge(initialTree.data, std::get<0>(link), std::get<1>(link));
        newTree.edges.push_back(newEdge);
        // add element to [(e', n')] vector
        updatedLinks.push_back(std::make_tuple(newEdge, std::get<1>(link)));
    }

    // synthesize new node from original node and [(e', n')] vector
    newTree.data = updateNode(initialTree.data, updatedLinks);

    return newTree;
}


// foldTree :: (l -> a -> a)
//          -> (n -> e -> a -> a)
//          -> a
//          -> RootedTree n e l
//          -> a
template<typename Node, typename Edge, typename Leaf, typename Accumulator>
Accumulator mlc_rooted_foldTree(
  std::function<Accumulator(Leaf,Accumulator)> foldLeaf,
  std::function<Accumulator(Node,Edge,Accumulator)> foldNode,
  Accumulator b,
  RootedTree<Node,Edge,Leaf> tree
){
    for (int i = 0; i < tree.children.size(); i++){
        auto child = tree.children[i];
        if (std::holds_alternative<RootedTree<Node, Edge, Leaf>>(child)) {
            const auto& child_tree = std::get<RootedTree<Node, Edge, Leaf>>(child);
            b = mlc_rooted_foldTree(foldLeaf, foldNode, b, child_tree);
            b = foldNode(child_tree.data, tree.edges[i], b);
        } else {
            const auto& child_leaf = std::get<Leaf>(child);
            b = foldLeaf(child_leaf, b);
        }
    }
    return b;
}

#endif
