#ifndef __ALGO_HPP__
#define __ALGO_HPP__

#include <mlc.hpp>


// Neighbor Joining Tree algorithm function
Tree<mlc::Unit, double, int> nj(const std::vector<std::string>& input) {
    // TODO: implement the function

    // returning a default tree for now
    return Tree<mlc::Unit, double, int>{mlc::unit, {}, {}};
}

// Multiple sequence alignment function
std::vector<std::string> align(const std::vector<std::string>& input) {
    // TODO: implement the function

    // returning the same input for now
    return input;
}

#endif
