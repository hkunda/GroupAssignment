#ifndef HUNGARIAN_H
#define HUNGARIAN_H

#include <vector>
#include <memory>

#include <Rcpp.h>
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::Dimension;

struct WeightedBipartiteEdge {
    const int left;
    const int right;
    const int cost;

    WeightedBipartiteEdge() : left(), right(), cost() {}
    WeightedBipartiteEdge(const int left, const int right, const int cost) : left(left), right(right), cost(cost) {}
};

// Given the number of nodes on each side of the bipartite graph and a list of edges, returns a minimum-weight perfect matching.
// If a matching is found, returns a length-n vector, giving the nodes on the right that the left nodes are matched to.
// If no matching exists, returns an empty vector.
// (Note: Edges with endpoints out of the range [0, n) are ignored.)
const std::vector<int> hungarianMinimumWeightPerfectMatching(int n, const std::vector<WeightedBipartiteEdge> allEdges);

/*
 * Used to compute the optimal assignment of students given a (partial) ranking of students.
 * Inputs:
 *  @rankings [NumericMatrix]: a matrix of rankings for each student. The ith row corresponds to student i,
 *  and the jth column corresponds to the student in position j on any student's ranking list. If there are
 *  N students, then each student should have an ID in the range [1,N], which is used both as the row index
 *  (to refer to a student's rankings) and entries of the matrix. For example, if rankings[i,j] = k, then
 *  student i has placed student k in position j on i's ranking.
 *
 *  @students [int]: the number of students. This should be equal to the number of rows in @rankings.
 *
 *  @leaders [int]: the number of leaders to pick; equivalently, the number of groups to form.
 *
 *  @minGroupSize [int]: the minimum number of students in a group. NOTE: this value does include the group
 *  leader. So if minGroupSize = 4, then every group must consist of 1 leader and at least 3 other students.
 *
 *  @maxGroupSize [int]: the maximum number of students in a group. NOTE: this value does include the group
 *  leader. So if maxGroupSize = 6, then every group must consist of 1 leader and at most 5 other students.
 *
 *  @missingCost [int]: the cost of assigninga student to a leader who isn't on their rankings.
 *  The above has not yet been implemented (issues with default value).
 *
 * Outputs:
 *  Returns a NumericVector of length students+1. The ith element of the output is the id of the student
 *  leading the group that student i is in (for 0 <= i <= N). The N+1st element of the output is the total
 *  cost of the optimal matching.
 */

// NumericVector optimalAssignment(NumericMatrix rankings, int students, int leaders, int minGroupSize, int maxGroupSize, int missingCost);
NumericVector optimalAssignment(NumericMatrix rankings, int students, int leaders, int minGroupSize, int maxGroupSize);

#endif // HUNGARIAN_H
