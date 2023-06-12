suppressMessages(require(ape))

randomTree <- function(n){
    x <- rtree(n)
    x$node.label <- paste0("n", 1:x$Nnode)
    x
}
