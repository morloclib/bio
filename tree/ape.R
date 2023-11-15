suppressMessages(require(ape))

random_rooted_tree <- function(n){
    x <- rtree(n, rooted = TRUE)
    x$node.label <- paste0("n", 1:x$Nnode)
    x
}

random_unrooted_tree <- function(n){
    x <- rtree(n, rooted = FALSE)
    x$node.label <- paste0("n", 1:x$Nnode)
    x
}
