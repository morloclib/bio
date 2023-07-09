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
#include <iostream>

#include <mlc.hpp>

template <typename Node, typename Edge, typename Leaf>
struct Tree {
    Node data;
    std::vector<std::variant<Tree<Node, Edge, Leaf>, Leaf>> children;
    std::vector<Edge> edges;
};

template <typename Node, typename Edge, typename Leaf>
int mlc_count_nodes(Tree<Node, Edge, Leaf> tree) {
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
        std::vector<Node> nodes,
        std::vector<Leaf> leafs,
        std::vector<std::vector<Edge>> edges,
        std::vector<std::vector<int>> child_indices) {
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
Tree<Node, Edge, Leaf> mlc_pack_tree(std::tuple<std::vector<Node>,
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
    return mlc_pack_tree_r(0, nodes, leafs, node_edges, child_indices);
}


template <typename Node, typename Edge, typename Leaf>
std::tuple<std::vector<Node>, std::vector<std::tuple<int, int, Edge>>, std::vector<Leaf>>
mlc_unpack_tree(const Tree<Node, Edge, Leaf> &tree)
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

// node :: Tree n e l -> n
template <typename Node, typename Edge, typename Leaf>
Node mlc_node(Tree<Node, Edge, Leaf> tree){
    return tree.data;
}

// subtrees :: Tree n e l -> [(e, Tree n e l)]
template <typename Node, typename Edge, typename Leaf>
std::vector<std::tuple<Edge, Tree<Node, Edge, Leaf>>> mlc_subtrees(Tree<Node, Edge, Leaf> tree){
    std::vector<std::tuple<Edge, Tree<Node, Edge, Leaf>>> subtrees;
    for(std::size_t i = 0; i < tree.children.size(); i++){
        if (std::holds_alternative<Tree<Node,Edge,Leaf>>(tree.children[i])) {
            subtrees.push_back(std::make_tuple(
              tree.edges[i],
              std::get<Tree<Node,Edge,Leaf>>(tree.children[i])
            ));
        }
    }
    return subtrees;
}

// childLeafs :: Tree n e l -> [l]
template <typename Node, typename Edge, typename Leaf>
std::vector<std::tuple<Edge, Leaf>> mlc_childLeafs(Tree<Node, Edge, Leaf> tree){
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
Tree<Node, Edge, NewLeaf> mlc_mapLeaf(
    std::function<NewLeaf(Leaf)> func, 
    Tree<Node, Edge, Leaf> tree
) {
    Tree<Node, Edge, NewLeaf> newTree;
    newTree.data = tree.data;
    newTree.edges = tree.edges;

    for (const auto& child : tree.children) {
        if (std::holds_alternative<Leaf>(child)) {
            // It's a leaf. Apply the function to it and add it to newTree.
            Leaf oldLeaf = std::get<Leaf>(child);
            NewLeaf newLeaf = func(oldLeaf);
            newTree.children.push_back(newLeaf);
        }
        else if (std::holds_alternative<Tree<Node, Edge, Leaf>>(child)) {
            // It's a tree. Recursively call mapLeafs on it and add it to newTree.
            const Tree<Node, Edge, Leaf>& oldSubtree = std::get<Tree<Node, Edge, Leaf>>(child);
            Tree<Node, Edge, NewLeaf> newSubtree = mlc_mapLeaf(func, oldSubtree);
            newTree.children.push_back(newSubtree);
        }
    }

    return newTree;
}


// treeBy :: ([b] -> Tree n e Int)
//        -> [(l, b)]
//        -> Tree n e l

// The list xs contains pairs of values. The first element in the pair is the
// value that will ultimately be stored in the tree's leafs. The second element
// is data that is used to generate the tree. The list xs is split into a list
// of B's that is passed to the buildTree function. The output tree contains
// integer indices as leaves. These indices are used to lookup the appropriate
// leaf in xs and replace the list in the tree.
template<typename Node, typename Edge, typename Leaf, typename B>
Tree<Node, Edge, Leaf> mlc_treeBy(
    std::function<Tree<Node, Edge, int>(std::vector<B>)> buildTree,
    std::vector<std::tuple<Leaf, B>> xs
) {
    // Extract the second element of each pair in xs to pass to buildTree.
    std::vector<B> buildTreeArgs;
    for (auto& pair : xs) {
        buildTreeArgs.push_back(std::get<1>(pair));
    }

    // Build the initial tree.
    Tree<Node, Edge, int> indexTree = buildTree(buildTreeArgs);

    Tree<Node, Edge, Leaf> finalTree;

    // STUB

    // // Replace indices with their corresponding leafs
    // Tree<Node, Edge, Leaf> finalTree = mlc_mapLeaf(
    //     [&](int index) { return std::get<0>(xs[index]); },
    //     indexTree
    // );

    return finalTree;
}


// push :: (a -> n -> n')
//      -> (a -> n' -> e -> n -> (e', n'))
//      -> (a -> n' -> e -> l -> (e', l'))
//      -> (Tree n e l -> a -> a)
//      -> a
//      -> Tree n e l
//      -> Tree n' e' l'
template<typename Accumulator, typename Node, typename Edge, typename Leaf, typename NodePrime, typename EdgePrime, typename LeafPrime>
Tree<NodePrime, EdgePrime, LeafPrime> mlc_push_val(
    std::function<NodePrime(Accumulator, Node)> handleRoot,
    std::function<std::tuple<EdgePrime, NodePrime>(Accumulator, NodePrime, Edge, Node)> alterChildNode,
    std::function<std::tuple<EdgePrime, LeafPrime>(Accumulator, NodePrime, Edge, Leaf)> alterLeaf,
    std::function<Accumulator(Tree<Node, Edge, Leaf>, Accumulator)> updateAccumulator,
    Accumulator a,
    Tree<Node, Edge, Leaf> initialTree
) {
    // Implementation goes here
    return Tree<NodePrime, EdgePrime, LeafPrime>{};
}


// push :: (n -> n')
//      -> (n' -> e -> n -> (e', n'))
//      -> (n' -> e -> l -> (e', l'))
//      -> Tree n e l
//      -> Tree n' e' l'
template<typename Node, typename Edge, typename Leaf, typename NodePrime, typename EdgePrime, typename LeafPrime>
Tree<NodePrime, EdgePrime, LeafPrime> mlc_push(
    std::function<NodePrime(Node)> handleRoot,
    std::function<std::tuple<EdgePrime, NodePrime>(NodePrime, Edge, Node)> alterChildNode,
    std::function<std::tuple<EdgePrime, LeafPrime>(NodePrime, Edge, Leaf)> alterLeaf,
    Tree<Node, Edge, Leaf> initialTree
) {
    // Implementation goes here
    return Tree<NodePrime, EdgePrime, LeafPrime>{};
}


// pull :: (l -> (a, n'))
//      -> (n -> e -> a -> n' -> e')
//      -> (n -> [(e', n', a)] -> (a, n'))
//      -> Tree n e l
//      -> Tree n' e' l
template<typename Leaf, typename A, typename NodePrime, typename Node, typename Edge, typename EdgePrime>
Tree<NodePrime, EdgePrime, Leaf> mlc_pull_val(
    std::function<std::tuple<A, NodePrime>(Leaf)> handleLeaf,
    std::function<EdgePrime(Node, Edge, A, NodePrime)> updateEdge,
    std::function<std::tuple<A, NodePrime>(Node, std::vector<std::tuple<EdgePrime, NodePrime, A>>)> updateNode,
    Tree<Node, Edge, Leaf> initialTree 
) {
    // Implementation goes here
    return Tree<NodePrime, EdgePrime, Leaf>{};
}


// template<typename Leaf, typename NodePrime, typename Node, typename Edge, typename EdgePrime>
// Tree<NodePrime, EdgePrime, Leaf> mlc_pull(
//     std::function<NodePrime(Leaf)> handleLeaf,
//     std::function<EdgePrime(Node, Edge, NodePrime)> updateEdge,
//     std::function<NodePrime(Node, std::vector<std::tuple<EdgePrime, NodePrime>>)> updateNode,
//     Tree<Node, Edge, Leaf> initialTree
// ) {
//     // Implementation goes here
//     return Tree<NodePrime, EdgePrime, Leaf>{};
// }

// pull :: (l -> n')
//      -> (n -> e -> n' -> e')
//      -> (n -> [(e', n')] -> n')
//      -> Tree n e l
//      -> Tree n' e' l
template<typename Leaf, typename NodePrime, typename Node, typename Edge, typename EdgePrime>
Tree<NodePrime, EdgePrime, Leaf> mlc_pull(
    std::function<NodePrime(Leaf)> handleLeaf,
    std::function<EdgePrime(Node, Edge, NodePrime)> updateEdge,
    std::function<NodePrime(Node, std::vector<std::tuple<EdgePrime, NodePrime>>)> updateNode,
    Tree<Node, Edge, Leaf> initialTree 
) {
    Tree<NodePrime, EdgePrime, Leaf> newTree;

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
        else if (std::holds_alternative<Tree<Node, Edge, Leaf>>(child)) {
            // It's a tree. Recursively call mlc_pull on it and add it to newTree.
            Tree<Node, Edge, Leaf> oldSubtree = std::get<Tree<Node, Edge, Leaf>>(child);
            Tree<NodePrime, EdgePrime, Leaf> newSubtree = mlc_pull(handleLeaf, updateEdge, updateNode, oldSubtree);
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


// foldTree :: (a -> a -> a)
//          -> (Tree n e l -> a)
//          -> a
//          -> Tree n e l
//          -> a
template<typename Accumulator, typename Node, typename Edge, typename Leaf>
Accumulator mlc_foldTree(
    std::function<Accumulator(Accumulator, Accumulator)> f,
    std::function<Accumulator(Tree<Node, Edge, Leaf>)> foldChild,
    Accumulator b,
    Tree<Node, Edge, Leaf> initialTree
) {
    // Implementation goes here
    return Accumulator{};
}

// alterTree :: (n -> [(e, Tree n e l)] -> [(e, Tree n e l)])
//           -> Tree n e l
//           -> Tree n e l
template<typename Node, typename Edge, typename Leaf>
Tree<Node, Edge, Leaf> mlc_alterTree(
    std::function<std::vector<std::tuple<Edge, Tree<Node, Edge, Leaf>>>(Node, std::vector<std::tuple<Edge, Tree<Node, Edge, Leaf>>>)> updateChildren,
    Tree<Node, Edge, Leaf> initialTree
) {
    // Implementation goes here
    return Tree<Node, Edge, Leaf>{};
}


#endif
