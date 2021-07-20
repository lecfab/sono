# Sorting nodes in networks

Here is a C++ tool to read a graph and permute its nodes according to various orders. The idea is that some algorithms may have different time and space efficiency depending on the ordering of nodes.


## Installation & compilation
`$ git clone https://github.com/lecfab/sono.git`

`$ make ord` to optimise the executable for production

`$ make ord DEBUG=1` to obtain a debugging executable (compile information, no optimisation, debug messages...)

## Running

`$ ./ord DATASET ORDER`

Type `$ ./ord --help` for more information.


#### Parameters
-   `DATASET`: a graph representation for nodes [0 to N-1]. It must be a text file where each line corresponds to a directed edge of the form `a b` (i.e. a SPACE b, with a and b long unsigned integers).

-   `ORDER`: currently available orders are optimise

    -   `original`:   the original order provided by the dataset (warning: it is usually not random)
    -   `rand`: random reordering of nodes
    -   `deg`:  sorted by decreasing total degree
    -   `deg+`, `deg-`: sorted by decreasing outgoing (resp. ingoing) degree
    -   `core`: sorted by degeneracy ordering (k-core pealing algorithm); note that for directed graphs, it corresponds to `core-` below
    -   `core+`, `core-`: k-core pealing using only out-edges (resp. in-edges)
    -   `icore`: inverted pealing, where nodes of higher degree are removed first
    -   `icore+`, `icore-`: inverted pealing, where nodes of higher out-degree (resp. in-degree) are removed first
    -   `dfs`, `bfs`: graph traversal starting from the node with smallest id
    -   `gorder`, `rcm`, `ldg`, modified `slashburn`, `minla`, `minloga` as described in [this paper](https://doi.org/10.5281/zenodo.4836230)
    -   `trianglesdpp`, `trianglesdpm`: orderings obtained by greedy optimisation of sum d+d+ (resp. d+d-) for each node, which is the approximate complexity of specific triangle-counting algorithms

#### Options
-   `-o FILE`: output file in which the order will be printed (the new ID of the node with the old ID u is on line u)
-   `-d`: flag to specify that the graph is directed (it <u>can</u> have both `a b` and `b a` edges)
-   `-u` (default): flag to specify that the graph is undirected (it <u>cannot</u> have both `a b` and `b a` edges). See below how to make a directed graph undirected.

#### Output
A file of N lines where line i contains the rank of node i according to ORDER.


## Tools
#### Undirected graphs
If your INPUT graph is directed and have opposite edges (eg. both `a b` and `b a` are in the edge list), we provide a tool to make it undirected. The OUTPUT edge list will only contain `a b`, while `b a` not written to reduce file size.
Here is the way to do so:

`$ make undirect`

`$ ./undirect INPUT OUTPUT`

#### Rank edges
The `ord` program outputs a rank: for N nodes, it contains N lines, where line i is the new rank of node i.

Here is a tool to renumber an edge list INPUT according to a rank RANK.

`$ make rankedges`

`$ ./rankedges INPUT RANK OUTPUT`

Note that RANK should be the output of `ord` applied on dataset INPUT.

## Contributors

Fabrice Lécuyer (fabrice.lecuyer@lip6.fr)

2021
