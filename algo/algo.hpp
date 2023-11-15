#ifndef __MORLOC_BIO_ALGO_HPP__
#define __MORLOC_BIO_ALGO_HPP__

#include "Eigen/Dense"
#include "mlccpptypes/prelude.hpp"
#include "mlccpptypes/bio.hpp"
#include <map>
#include "math.h"
#include <iostream>
using namespace std;

// countKmers :: k:Int -> s:Str -> Map Str Int
std::map<std::string,int> countKmers(int k, std::string seq){
    std::map<std::string,int> kmers;
   
    if(k <= 0){
        throw std::invalid_argument("k must be an integer greater than or equal to 1");
    }

    for(int i = 0; (i + k) <= seq.size(); i++){
        std::string kmer = seq.substr(i, k);
        if(kmers.find(kmer) != kmers.end()){
            kmers[kmer]++;
        } else {
            kmers[kmer] = 1;
        }
    }

    return kmers;
}

// kmerDistance :: Map Str Int -> Map Str Int -> Real
//
// distance = sqrt(sum( (x_i - y_i)^2 / (x_i + y_i) ))
double kmerDistance(const std::map<std::string,int>& x, std::map<std::string,int> y){
    double square_distance = 0.0;

    // Iterate through kmers in sequence x, find distance to y
    for (const auto& pair : x){
        std::string key = pair.first;
        int xcount = pair.second;
        int ycount = 0; 
        if(y.find(key) != y.end()){
            ycount = y[key]; 
        }
        square_distance += (xcount - ycount) * (xcount - ycount) / (xcount + ycount);
    }

    // Iterate through kmers in y, if they are missing from x, then they would
    // not have been accounted for in the prior loop, and the distance is simply
    // y^2
    for (const auto& pair : y){
        std::string key = pair.first;
        int ycount = pair.second;
        if(x.find(key) == x.end()){
            square_distance += ycount * ycount;
        }
    }

    return std::sqrt(square_distance);
}

/*

 1    2  3    4  5
(A, ((B, C), (D, E)))

6 (2, 3)

 1   6    4  5
(A, (BC, (D, E)))

7 (4, 5)

 1   6    7
(A, (BC, (DE)))

8 (6, 7)

 1  8
(A, BCDE)

0 (1, 8)

*/



typedef struct IndexPair{
  size_t a;
  size_t b;
  double dist;
} IndexPair;

RootedTree<mlc::Unit, double, int> edge_map_to_tree(std::vector<IndexPair> nodes){
    std::vector<RootedTree<mlc::Unit, double, int>> subtrees;
    std::vector<int> indices(nodes.size() + 1, -1);
    int Nnodes = 0;
    for (const auto& node : nodes) {
        int a = node.a;
        int b = node.b;
        if (indices[a] == -1) {
            subtrees.emplace_back(RootedTree<mlc::Unit, double, int>{});
            subtrees.back().data = mlc::Unit{};
            subtrees.back().children.emplace_back(a);
            indices[a] = subtrees.size() - 1;
        }
        if (indices[b] == -1) {
            subtrees.emplace_back(RootedTree<mlc::Unit, double, int>{});
            subtrees.back().data = mlc::Unit{};
            subtrees.back().children.emplace_back(b);
            indices[b] = subtrees.size() - 1;
        }
        RootedTree<mlc::Unit, double, int> new_tree{};
        new_tree.data = mlc::Unit{};
        new_tree.children.emplace_back(subtrees[indices[a]]);
        new_tree.children.emplace_back(subtrees[indices[b]]);
        new_tree.edges.emplace_back(node.dist);
        subtrees.emplace_back(new_tree);
        indices[a] = indices[b] = subtrees.size() - 1;
        Nnodes++;
    }
    return subtrees.back();
}


// upgmaFromDist :: Matrix Real -> RootedTree () Real Int This is a naive cubic
// time algorithm. Quadratic time algorithms are possible, this implementation
// just serves as a baseline.
RootedTree<mlc::Unit, double, int> upgmaFromDist(
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat
){
    int Nleafs = mat.cols();
    int Nnodes = 0;
    std::vector<IndexPair> nodes;
    std::vector<int> indices(Nleafs);
    for(int i = 0; i < indices.size(); i++){
      indices[i] = i;
    }

    bool not_done = true;
    int loop_count = 0;
    while(not_done){
      std::cerr << "loop " << loop_count << std::endl; loop_count++;
      double min_dist = 0;
      int min_i = 0;
      int min_j = 0;
      // Find closest pair of taxa
      for(int i = 0; i < Nleafs; i++){
        if(indices[i] < 0){
          for(int j = 0; j < i; j++){
            if(indices[j] < 0){
              not_done = false;
              double ij_dist = mat(i,j);
              if(ij_dist > -1 * min_dist){
                min_dist = ij_dist;
                min_i = i;
                min_j = j;
              }
            }
          }
        }
      }
      if(not_done){
        indices[min_i] = -1; // this row will no longer be used
        Nnodes++;
        indices[min_j] = Nnodes + Nleafs;

        IndexPair newNode;
        newNode.a = min_i;
        newNode.b = min_j;
        newNode.dist = min_dist;
        nodes.push_back(newNode);

        //      j
        //      |
        //    01234567             01834567
        //   1 0                  1 0
        //   2 -0                 8 +0
        //   3  |0                3  +*
        // i-4**m*0     [1 / 9]   4*****
        //   5  |* 0              5  +* 0
        //   6  |*  0             6  +*  0
        //   7  |*   0            7  +*   0
        int new_node = min_j; // to save space, we will reuse the min_j'th col
        for(int i = 0; i < min_j; i++){
          if(indices[i] >= 0){
            mat(new_node,i) = mat(min_j,i) + mat(min_i,i);
          }
        }
        for(int i = min_j + 1; i < Nleafs; i++){
          if(indices[i] >= 0){
            mat(i,new_node) = mat(i,min_j) + mat(min_i,i);
          }
        }
      }
    }
    std::cerr << "loop complete" << std::endl;

    RootedTree<mlc::Unit, double, int> finalTree = edge_map_to_tree(nodes);
    return finalTree;
}

#endif
