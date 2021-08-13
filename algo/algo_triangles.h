#ifndef ALGO_TRIANGLES_H
#define ALGO_TRIANGLES_H

#include "../utils/tools.h"
struct Edgelist;
struct Adjlist;
struct Badjlist;

struct pair_hash
{
  template <class T1, class T2>
  std::size_t operator() (const std::pair<T1, T2> &pair) const {
    auto h1 = std::hash<T1>()(pair.first);
    auto h2 = std::hash<T2>()(pair.second);
    h1 ^= h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);
    return h1;
  }
};

void triangle_complexities(const Badjlist &g);

ull count_triangles_hash(const Badjlist &g);
ull count_triangles_bool(const Badjlist &g, bool uOut_iter, bool vOut_iter);
ull count_triangles_dpp(const Badjlist &g);
ull count_triangles_dmm(const Badjlist &g);
ull count_triangles_dpm(const Badjlist &g);
ull count_triangles_dpm_indep(const Badjlist &g);
ull count_triangles_dpm_indep_vectint(const Badjlist &g);
ull count_triangles_dichotomy(const Badjlist &g);
bool is_neighOut_dichotomy(const Badjlist &g, const ul &u, const ul &v, ull &c);



double burden_equal(const Adjlist &g);
ul burden_equal_int(const Adjlist &g);
ul burden_loop(const Adjlist &g, std::vector<ul> &degrees, bool permutation=false);
ul burden_periphery(const Adjlist &g);
ul burden_periphery(const Badjlist &g);
ul burden_permutation(const Adjlist &g);
ul burden_permutation(const Badjlist &g);

ull count_cliques(const Badjlist &g, ul k);
ull count_cliques_parallel(const Badjlist &g, ul k, int p);
ull count_cliques_parallel_(Badjlist &g, ul k, int p);
ull count_cliques_parallel_edges(const Badjlist &g, const Edgelist &h, ul k, int p);
ull count_cliques_5(Badjlist &g);
ull count_cliques_5_(Badjlist &g);

#endif
