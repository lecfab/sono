#ifndef ORDER_TRIANGLEDPP_H
#define ORDER_TRIANGLEDPP_H


#include "../utils/tools.h"

struct Adjlist;
struct ContinuousRank;

// objective functions
ul sum_square_dpp(const Adjlist &g, const std::vector<ul> &vec);
ul sum_square_dpm(const Adjlist &g, const std::vector<ul> &vec);
ul sum_dpm2(const Adjlist &g, const std::vector<ul> &vec);
ul sum_dpm_ab(const Adjlist &g, const std::vector<ul> &vec, double a, double b);
ul sum_mindpp(const Adjlist &g, const std::vector<ul> &vec);

// search for position in neighbourhood
std::pair<ul, long long int> best_position_dpp(const Adjlist &g, const ul &u, const std::vector<ul> &neighbours, std::vector<ul> &dp, ul i0);
std::pair<ul, long long int> best_position_dpm(const Adjlist &g, const ul &u, const std::vector<ul> &neighbours, std::vector<ul> &dp, ul i0);
std::pair<ul, long long int> best_position_dpm2(const Adjlist &g, const ul &u, const std::vector<ul> &neighbours, std::vector<ul> &dp, ul i0);
std::pair<ul, long long int> best_position_dpm_ab(const Adjlist &g, const ul &u, const std::vector<ul> &neighbours, std::vector<ul> &dp, ul i0, double a, double b);
std::pair<ul, long long int> best_position_mindpp(const Adjlist &g, const ul &u, const std::vector<ul> &neighbours, std::vector<ul> &dp, ul i0);

// main optimisations
std::vector<ul> place_neighbour(const Adjlist &g, decltype(sum_square_dpm) &scoring, decltype(best_position_dpm) &best_position);
std::vector<ul> place_neighbour_dpp(const Adjlist &g);
std::vector<ul> place_neighbour_dpm(const Adjlist &g);
std::vector<ul> place_neighbour_dpm2(const Adjlist &g);
std::vector<ul> place_neighbour_mindpp(const Adjlist &g);
// std::vector<ul> optimise_edges(const Adjlist &g);

// compute current out-going degree
std::vector<ul> compute_dp(const Adjlist &g, const ContinuousRank &cRank);

std::list<ul>::iterator insert_into_chain(std::list<ul> &chain, std::list<ul>::iterator pointer_u, std::list<ul>::iterator pointer_next, std::vector<double> &altitude);

std::vector<ul> rank_from_chain(const std::list<ul> &chain);

#endif
