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

// node :: Tree n e l -> n
template <typename Node, typename Edge, typename Leaf>
Node mlc_node(Tree<Node, Edge, Leaf> tree){
    return tree.data;
}

// kids :: Tree n e l -> [(e, Tree n e l)]
template <typename Node, typename Edge, typename Leaf>
std::vector<std::tuple<Edge, std::variant<Tree<Node, Edge, Leaf>, Leaf>>> mlc_kids(Tree<Node, Edge, Leaf> tree){
    return tree.children;
}

// childLeafs :: Tree n e l -> [l]
template <typename Node, typename Edge, typename Leaf>
std::vector<Leaf> mlc_childLeafs(Tree<Node, Edge, Leaf> tree){
    std::vector<Leaf> leafs;
    for (const auto& child : tree.children) {
        if (std::holds_alternative<Leaf>(child)) {
            leafs.push_back(child);
        }
    }
    return leafs;
}

template <typename Node, typename Edge, typename Leaf, typename NewLeaf>
Tree<Node, Edge, NewLeaf> mapLeafs(
    std::function<NewLeaf(Leaf)> func, 
    const Tree<Node, Edge, Leaf>& tree
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
            Tree<Node, Edge, NewLeaf> newSubtree = mapLeafs(func, oldSubtree);
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
Tree<Node, Edge, Leaf> treeBy(
    std::function<Tree<Node, Edge, int>(std::vector<B>)> buildTree,
    std::vector<std::pair<Leaf, B>> xs
) {
    // Extract the second element of each pair in xs to pass to buildTree.
    std::vector<B> buildTreeArgs;
    for (auto& pair : xs) {
        buildTreeArgs.push_back(pair.second);
    }

    // Build the initial tree.
    Tree<Node, Edge, int> indexTree = buildTree(buildTreeArgs);

    // Replace each leaf index with the corresponding leaf value from xs.
    Tree<Node, Edge, Leaf> finalTree;
    // Assuming the Tree structure has a function `traverse` which traverse through the tree and replace leaf values
    indexTree.traverse([&](int index) {
        return xs[index].first;
    }, finalTree);

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
Tree<NodePrime, EdgePrime, LeafPrime> push(
    std::function<NodePrime(Accumulator, Node)> handleRoot,
    std::function<std::pair<EdgePrime, NodePrime>(Accumulator, NodePrime, Edge, Node)> alterChildNode,
    std::function<std::pair<EdgePrime, LeafPrime>(Accumulator, NodePrime, Edge, Leaf)> alterLeaf,
    std::function<Accumulator(Tree<Node, Edge, Leaf>, Accumulator)> updateAccumulator,
    Accumulator a,
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
Tree<NodePrime, EdgePrime, Leaf> pull(
    std::function<std::pair<A, NodePrime>(Leaf)> handleLeaf,
    std::function<EdgePrime(Node, Edge, A, NodePrime)> updateEdge,
    std::function<std::pair<A, NodePrime>(Node, std::vector<std::tuple<EdgePrime, NodePrime, A>>)> updateNode,
    Tree<Node, Edge, Leaf> initialTree 
) {
    // Implementation goes here
    return Tree<NodePrime, EdgePrime, Leaf>{};
}

// foldTree :: (a -> a -> a)
//          -> (Tree n e l -> a)
//          -> a
//          -> Tree n e l
//          -> a
template<typename Accumulator, typename Node, typename Edge, typename Leaf>
Accumulator foldTree(
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
Tree<Node, Edge, Leaf> alterTree(
    std::function<std::vector<std::pair<Edge, Tree<Node, Edge, Leaf>>>(Node, std::vector<std::pair<Edge, Tree<Node, Edge, Leaf>>>)> updateChildren,
    Tree<Node, Edge, Leaf> initialTree
) {
    // Implementation goes here
    return Tree<Node, Edge, Leaf>{};
}


#endif
