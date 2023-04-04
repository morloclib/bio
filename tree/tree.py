import dendropy

#  tree.mlc_pack_tree_dendropy(([1,2,3],[(0,1,"hi"),(0,2,"bi")],["sussy", "lucy"]))

#  packTree Py   :: pack   => ([n], [("int", "int", e)], [l]) -> "dendropy" n e l
def mlc_pack_tree_dendropy(x):
    nodes = x[0]
    leafs = x[2]
    edge_map = [[] for _ in range(len(x[1]))]
    for (p, c, _) in x[1]:
        edge_map[p].append(c)

    edges = [e for (_, _, e) in x[1]]

    taxon_namespace = dendropy.TaxonNamespace(leafs)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)

    return pack_tree_dendropy_r(tree, nodes, leafs, edges, edge_map, 0)

def pack_tree_dendropy_r(tree, nodes, leafs, edges, edge_map, index):
    if(index >= len(nodes)):
        # this is a leaf
        tree.label = leafs[index - len(nodes)]
    else:
        # this is a node
        tree.label = nodes[index]
        children = []
        for child_index in edge_map[index]:
            if child_index >= len(nodes):
                child_edge = edges[child_index - len(nodes)]
            else:
                child_edge = edges[child_index]
            child_node = dendropy.Node(edge_length=child_edge)
            child_tree = pack_tree_dendropy_r(child_node, nodes, leafs, edges, edge_map, child_index)
            children.append(child_tree)
        tree.set_child_nodes(children)

    return tree


#  unpackTree Py :: unpack => "dendropy" n e l -> ([n], [("int", "int", e)], [l])
def mlc_unpack_tree_dendropy(tree):
    pass

#  taxon_namespace = dendropy.TaxonNamespace(["A", "B", "C", "D",])
#  tree = dendropy.Tree(taxon_namespace=taxon_namespace)
#
#
#  # Create and add a new child node to the seed node,
#  # assigning it an edge length:
#  #
#  #     (seed)
#  #      /
#  #     /
#  #    ch1
#  #
#  ch1 = tree.seed_node.new_child()
#  ch1.edge.length = 1
#
#  # Can also assign edge length on construction:
#  #
#  #     (seed)
#  #      / \
#  #     /   \
#  #   ch1   ch2
#  #
#  ch2 = tree.seed_node.new_child(edge_length=1)
#
#  # Can also add an existing node as child
#  #
#  #       (seed)
#  #       /   \
#  #      /     \
#  #    ch1     ch2
#  #   /  \     /  \
#  #  ch3 ch4  ch5 ch6
#  ch3 = dendropy.Node(edge_length=1)
#  ch4 = dendropy.Node(edge_length=2)
#  ch1.add_child(ch3)
#  ch1.add_child(ch4)
#  ch5 = dendropy.Node(edge_length=1)
#  ch6 = dendropy.Node(edge_length=2)
#  # Note: this clears/deletes existing child nodes before adding the new ones;
#  ch2.set_child_nodes([ch5, ch6])
#
#  # Assign taxa
#  ch3.taxon = taxon_namespace.get_taxon("A")
#  ch4.taxon = taxon_namespace.get_taxon("B")
#  ch5.taxon = taxon_namespace.get_taxon("C")
#  ch6.taxon = taxon_namespace.get_taxon("D")
#
#  print(tree.as_string("newick"))
#  print(tree.as_ascii_plot())
