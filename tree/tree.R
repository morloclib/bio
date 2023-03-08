# packTree R   :: pack   => ([n], [("int", "int", e)], [l]) -> "phylo" n e l
mlc_pack_tree <- function(x){
    nodes <- x[[1]]
    edges <- lapply(x[[2]], function(e) e[[3]])
    leafs <- x[[3]]

    edge_map <- lapply(x[[2]], function(e) { c(e[[1]], e[[2]])})
    edge_map[,1] <- edge_map[,1] + length(leafs) + 1
    edge_map[,2] <- iflese(edge_map[,2] < length(nodes), edge_map[,2] + length(leafs) + 1, edge_map[,2] - length(nodes) + 1)

    tree <- list()
    class(tree) <- "phylo"
    tree$edge <- edge_map
    tree$tip.label <- leafs
    tree$node.label <- nodes
    tree$edge.length <- edges
    tree
}

# unpackTree R :: unpack => "phylo" n e l -> ([n], [("int", "int", e)], [l])
mlc_unpack_tree <- function(tree){
    if(is.null(tree$edge.length)){
        tree$edge.length <- lapply(1:nrow(tree$edge), function(i) "")
    }
    if(is.null(tree$node.label)){
        tree$node.label <- lapply(1:nrow(tree$edge), function(i) "")
    }
    nodes <- tree$node.label
    leafs <- tree$tip.label
    edges <- tree$edge.length
    N <- length(tree$tip.label)

    # morloc edge maps have 0 as root, nodes as indices (0..N-1),
    # and leafs as indices (N..(L+N-1)).
    # where
    #   N is number of nodes
    #   L is number of leafs
    # R phylo objects are the opposite, so here I need to reverse
    edge_map <- lapply(1:nrow(tree$edge), function(i) {
        node_idx <- tree$edge[i, 1] - N - 1
        leaf_idx <- if(tree$edge[i, 2] > N){
            tree$edge[i,2] - N - 1
        } else {
            tree$edge[i,2] + tree$Nnode - 1
        }
        list(node_idx, leaf_idx, edges[i])
    })

    list(nodes, edge_map, leafs)
}
