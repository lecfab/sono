// See header for general documentation

#include <functional>
#include <unordered_set>
#include <set>
#include "algo_triangles.h"
#include "../utils/adjlist.h"
#include "../utils/edgelist.h"
#include <iomanip>
#include <omp.h>

using namespace std;


double burden_equal(const Adjlist &g) {
  double m = g.e / g.edge_factor;
  return m*m / g.n; // = (m/n)² * n
}

ul burden_equal_int(const Adjlist &g) {
  ul m = g.e / g.edge_factor;
  ul burden = m / g.n, nodes_above = m % g.n;
  return burden*burden * (g.n - nodes_above) + (burden+1)*(burden+1) * nodes_above;
}

ul burden_loop(const Adjlist &g, vector<ul> &degrees, bool permutation) {
  sort(degrees.begin(), degrees.end());
  ul dpp = 0, n = g.n, m = g.e / g.edge_factor;
  for(ul i=0; i<g.n; ++i) {
    ul burden = min(degrees[i], m / n); // the node takes as much as it can but no more than the equal remaining burden
    burden = min(burden, i + n*(!permutation)); // in a permutation, the ith-last node can only have i out-neighbours
    dpp += burden * burden;
    m -= burden;
    n --;
  }
  return dpp;
}


ul burden_periphery(const Adjlist &g) {
  vector<ul> degrees; degrees.reserve(g.n);
  for(ul u=0; u<g.n; ++u) degrees.push_back(g.get_degree(u));
  return burden_loop(g, degrees, false);
}
ul burden_periphery(const Badjlist &g) {
  vector<ul> degrees; degrees.reserve(g.n);
  for(ul u=0; u<g.n; ++u) degrees.push_back(g.get_deg(u)); // in+out degree
  return burden_loop(g, degrees, false);
}

ul burden_permutation(const Adjlist &g) {
  vector<ul> degrees; degrees.reserve(g.n);
  for(ul u=0; u<g.n; ++u) degrees.push_back(g.get_degree(u));
  return burden_loop(g, degrees, true);
}
ul burden_permutation(const Badjlist &g) {
  vector<ul> degrees; degrees.reserve(g.n);
  for(ul u=0; u<g.n; ++u) degrees.push_back(g.get_deg(u));
  return burden_loop(g, degrees, true);
}




void triangle_complexities(const Badjlist &g) {
  double m = g.e / g.edge_factor;
  ull minpp = 0, minpm = 0/*, minmp = 0*/, minmm = 0;
  for(ul u=0; u<g.n; ++u)
    for(auto &v : g.neighOut_iter(u)) {
      minpp += min(g.get_degOut(u), g.get_degOut(v));
      minpm += min(g.get_degOut(u), g.get_degIn(v));
      // minmp += min(g.get_degIn(u), g.get_degOut(v)); // impossible: uvw would be a cycle
      minmm += min(g.get_degIn(u), g.get_degIn(v));
    }
  Info("O(sum min(du+, dv+)) : "<< minpp << " \t| per edge: " << ((double) minpp)/m)
  Info("O(sum min(du+, dv-)) : "<< minpm << " \t| per edge: " << ((double) minpm)/m)
  // Info("O(sum min(du-, dv+)) : "<< minmp << " \t| per edge: " << ((double) minmp)/m)
  Info("O(sum min(du-, dv-)) : "<< minmm << " \t| per edge: " << ((double) minmm)/m)

  ull dpp = 0, dpm = 0, dmm = 0, dep = 0, dpm2 = 0;
  for(ul u=0; u<g.n; ++u) {
    dpp += g.get_degOut(u) * g.get_degOut(u);
    dpm += g.get_degOut(u) * g.get_degIn(u);
    dpm2 += pow(g.get_degOut(u) * g.get_degIn(u), 2);
    // if(g.get_degOut(u) * g.get_degIn(u)) { // debug for bipartite graph
    //   cout << u << " "<<g.get_degOut(u) <<" "<< g.get_degIn(u)<<endl;
    //   cout << "OUT:\t";
    //   for(auto &v: g.neighOut_iter(u)) cout << v << "\t";
    //   cout << endl << "IN:\t";
    //   for(auto &v: g.neighIn_iter(u)) cout << v << "\t";
    //   cout << endl;
    // }
    dmm += g.get_degIn(u) * g.get_degIn(u);
    dep += g.get_deg(u) * g.get_degOut(u);
  }

  Info("O(sum d+²) ≥ permut "<< burden_permutation(g)/m <<" ≥ periph "<< burden_periphery(g)/m <<" ≥ int " << burden_equal_int(g)/m <<" ≥ eq "<< burden_equal(g)/m)
  Info("O(sum d+²) : "<< dpp << " \t| per edge: " << ((double) dpp)/m)
  Info("O(sum d-²) : "<< dmm << " \t| per edge: " << ((double) dmm)/m)
  Info("O(sum d+d-) : "<< dpm << " \t| per edge: " << ((double) dpm)/m)
  Info("O(sum d*d+) : "<< dep << " \t| per edge: " << ((double) dep)/m)
  Info("O(sum d+²d-²) : "<< dpm2 << " \t| per edge: " << ((double) dpm2)/m)
}

ull count_triangles_hash(const Badjlist &g) {
  Info("Counting triangles with hash complexity O(sum min(du+, dv+))")
  double m = g.e / g.edge_factor;
  unordered_set<edge,pair_hash> edges; edges.reserve(g.e);
  // set<edge> edges;
  // unordered_set<ull> edges;
  for(ul u=0; u<g.n; ++u) {
    // ull offset_u = g.n*u;
    for(auto &v : g.neighOut_iter(u))
      edges.insert({u, v});
      // edges.emplace(edge(u, v));
      // edges.insert(offset_u + v);
  }
  Info("Fillled")
  ull t = 0, c = 0;
  for(ul u=0; u<g.n; ++u) {
    for(auto &v : g.neighOut_iter(u)) {
      if(g.get_degOut(v) < g.get_degOut(u)) {
        // ull offset_u = g.n*u;
        for(auto &w : g.neighOut_iter(v)) {
          t += edges.count({u, w});
          // t += edges.count(edge(u, w));
          // t += edges.count(offset_u + w);
          c++;
        }
      }
      else {
        // ull offset_v = g.n*v;
        for(auto &w : g.neighOut_iter(u)) {
          t += edges.count({v, w});
          // t += edges.count(edge(v, w));
          // t += edges.count(offset_v + w);
          c++;
        }
      }
    }
  }
  Info("Found triangles: "<< t << " \t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  return t;
}



inline bool is_neighOut_dichotomy(const Badjlist &g, const ul &u, const ul &v, ull &c) {
  Debug("start "<<u<<" "<<v)
  // if(g.get_degOut(u) < 8) { // linear pass for small lists // useless because we only use dichotomy for big neighbourhoods
  //   for(auto &w : g.neighOut_iter(u)) if(w == v) return true;
  //   return false;
  // } // dichotomy for bigger lists
  ul x = g.cd[u], z = g.cd[u+1];
  while(x < z) {
    c++;
    ul y = (z+x) / 2; Debug("dic "<< x<<" "<<y <<" "<< z)
    if     (g.adj[y] > v) z = y;
    else if(g.adj[y] < v) x = y+1;
    else return true;
  } Debug("false")
  return false;
}

ull count_triangles_dichotomy(const Badjlist &g) {
  // if(!g.are_neighbours_sorted()) {
  //   Alert("Adjacency lists are not sorted")
  //   exit(1);
  // }
  double m = g.e / g.edge_factor;
  ull t = 0, c = 0, c2 = 0, left=0, right=0;

  vector<ul> dicho_speedup; dicho_speedup.reserve(g.n); // estimation of dichotomy complexity versus linear search
  for(ul u=0; u<g.n; ++u) dicho_speedup.push_back( g.get_degOut(u) / (1. + 2.*log2(1+g.get_degOut(u))) );

  vector<bool> is_neighOut(g.n, false);
  for(ul u=0; u<g.n; ++u) {
    for(auto &v : g.neighOut_iter(u)) is_neighOut[v] = true; // set array
    for(auto &v : g.neighOut_iter(u)) {
      if(dicho_speedup[v] < g.get_degOut(u)) {
        left++;
        for(auto &w : g.neighOut_iter(v)) {
          c++; c2++;
          t += is_neighOut[w];
        }
      }
      else {
        Debug(g.get_degOut(v) << " donne " << g.get_degOut(u) << " ≤ " << dicho_speedup[v])
        right++;
        for(auto &w : g.neighOut_iter(u)) {
          c++;
          t += is_neighOut_dichotomy(g, v, w, c2);
        }
      }
    }
    for(auto &v : g.neighOut_iter(u)) is_neighOut[v] = false; // reset array
  }
  Info("Found triangles: "<< t << "\t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  Info("Complexit2: "<< c2 << " \t| per edge: " << ((double) c2) / m)
  Info("repartition: "<< left << " \t|    " << right)
  return t;
}

ull count_triangles_bool(const Badjlist &g, bool uOut_iter, bool vOut_iter) {
  double m = g.e / g.edge_factor;
  function<NeighIter(ul)> Out_iter = [&g] (ul u) -> NeighIter { return g.neighOut_iter(u); };
  function<NeighIter(ul)> In_iter = [&g] (ul u) -> NeighIter { return g.neighIn_iter(u); };
  auto u_iter = uOut_iter ? Out_iter : In_iter;
  auto v_iter = vOut_iter ? Out_iter : In_iter;

  ull t = 0, c = 0;
  vector<bool> is_neighOut(g.n, false);
  for(ul u=0; u<g.n; ++u) {
    for(auto &v : g.neighOut_iter(u)) is_neighOut[v] = true; // set array

    for(auto &v : u_iter(u)) {
      for(auto &w : v_iter(v)) {
        c ++;
        t += is_neighOut[w];
      }
    }

    for(auto &v : g.neighOut_iter(u)) is_neighOut[v] = false; // reset array
  }
  Info("Found triangles: "<< t << " \t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  return t;
}

ull count_triangles_dpm_indep(const Badjlist &g) {
  Info("Counting triangles whith complexity O(sum d+d-) -- Independant function")
  double m = g.e / g.edge_factor;

  ull t = 0, c = 0;
  vector<bool> is_neighOut(g.n, false);
  for(ul u=0; u<g.n; ++u) {
    for(auto &v : g.neighOut_iter(u)) is_neighOut[v] = true; // set array

    for(auto &v : g.neighOut_iter(u)) {
      for(auto &w : g.neighOut_iter(v)) {
        c ++;
        t += is_neighOut[w];
      }
    }

    for(auto &v : g.neighOut_iter(u)) is_neighOut[v] = false; // reset array
  }
  Info("Found triangles: "<< t << " \t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  return t;
}

ull count_triangles_dpm_indep_vectint(const Badjlist &g) {
  Info("Counting triangles whith complexity O(sum d+d-) -- Independant function")
  double m = g.e / g.edge_factor;

  ull t = 0, c = 0;
  vector<int> is_neighOut(g.n, false);
  for(ul u=0; u<g.n; ++u) {
    for(auto &v : g.neighOut_iter(u)) is_neighOut[v] = true; // set array

    for(auto &v : g.neighOut_iter(u)) {
      for(auto &w : g.neighOut_iter(v)) {
        c ++;
        t += is_neighOut[w];
      }
    }

    for(auto &v : g.neighOut_iter(u)) is_neighOut[v] = false; // reset array
  }
  Info("Found triangles: "<< t << " \t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  return t;
}

ull count_triangles_dpm(const Badjlist &g) {
    Info("Counting triangles whith complexity O(sum d+d-)")
    return count_triangles_bool(g, true, true);
  }

ull count_triangles_dpp(const Badjlist &g) {
  Info("Counting triangles whith complexity O(sum d+²)")
  return count_triangles_bool(g, false, true);
}

ull count_triangles_dmm(const Badjlist &g) {
  Info("Counting triangles whith complexity O(sum d-²)")
  return count_triangles_bool(g, true, false);
}

// ull count_cliques(const Badjlist &g, ul k) {
//   Info("Counting cliques for k="<< k)
//   double m = g.e / g.edge_factor;
//
//   ull t = 0, c = 0;
//   vector<ul> nodes; nodes.reserve(g.n);
//   for(ul u=0; u<g.n; ++u) nodes.push_back(u);
//   nodes.push_back(g.n + 2); // useless node added
//   vector<ul> parents(g.n, 0);
//
//   vector<const ul*> iters_end; iters_end.reserve(k);
//   vector<const ul*> iters; iters.reserve(k);
//   iters.push_back(&nodes[0]);
//   iters_end.push_back(&nodes.back());
//
//   while(true) {
//     // for(auto &uu: iters) cout << *uu <<"\t"; cout << endl;
//
//     if(iters.back() == iters_end.back()) {
//       Debug("No more neighbours for floor "<< iters.size())
//       if(iters.size() == 1) break;
//       iters.pop_back();
//       iters_end.pop_back();
//       for(auto &w : g.neighOut_iter(*iters.back())) parents[w] --;
//       ++iters.back();
//     }
//     else if(iters.size() < k-1) {
//       ul v  = *iters.back();
//       // Debug("New floor above "<< iters.size() << " based on "<<v)
//       if(parents[v] < iters.size()-1) { // not all previous nodes are parents
//         Debug("Unsatisfying: "<<parents[v])
//         ++iters.back();
//         continue;
//       }
//       for(auto &w : g.neighOut_iter(v)) parents[w] ++;
//       iters.push_back(g.neighOut_beg(v));
//       iters_end.push_back(g.neighOut_end(v));
//     }
//     else {
//       ul x = *iters.back();
//       ++iters.back();
//       if(parents[x] < iters.size()-1) continue;
//       for(auto &y : g.neighOut_iter(x)) {
//         c++;
//         if(parents[y] == iters.size()-1) {
//           t++; // is neighbour of u, v, w (and x)
//           // Info(*iters[0]<<" "<<*iters[1]<<" "<<*iters[2]<<" "<<*iters[3]<<" "<<y)
//           // if(t > 1000000000) return t;
//         }
//       }
//     }
//   }
//   Info("Found "<<k<<"-cliques: "<< t << " \t| per edge: " << ((double) t) / m)
//   Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
//   return t;
// }
ull count_cliques(const Badjlist &g, ul k) { return count_cliques_parallel(g, k, 1); }
ull count_cliques_parallel(const Badjlist &g, ul k, int p) {
  omp_set_num_threads(p);
  Info("Counting cliques for k="<< k<<" with threads "<<p)
  double m = g.e / g.edge_factor;

  ull t = 0, c = 0;

  #pragma omp parallel reduction(+ : t, c)
  {
    vector<ul> parents(g.n, 0);
    #pragma omp for schedule(dynamic, 1)
    for(ul u=0; u<g.n; ++u) {
      for(auto &w : g.neighOut_iter(u))
        parents[w] ++;
      vector<const ul*> iters_end; iters_end.reserve(k-1);
      vector<const ul*> iters; iters.reserve(k-1);
      iters.push_back(g.neighOut_beg(u));
      iters_end.push_back(g.neighOut_end(u));
      while(true) {
        if(iters.back() == iters_end.back()) {
          // Debug("No more neighbours for floor "<< iters.size())
          if(iters.size() == 1) break;
          iters.pop_back();
          iters_end.pop_back();
          for(auto &w : g.neighOut_iter(*iters.back())) parents[w] --;
          ++iters.back();
        }
        else if(iters.size() < k-2) {
          ul v  = *iters.back();
          if(parents[v] < iters.size()) { // not all previous nodes are parents
            // Debug("Unsatisfying: "<<parents[v])
            ++iters.back();
            continue;
          }
          for(auto &w : g.neighOut_iter(v)) parents[w] ++;
          iters.push_back(g.neighOut_beg(v));
          iters_end.push_back(g.neighOut_end(v));
        }
        else {
          ul x = *iters.back();
          ++iters.back();
          if(parents[x] < iters.size()) continue;
          for(auto &y : g.neighOut_iter(x)) {
            c++;
            if(parents[y] == iters.size()) {
              t++; // is neighbour of u, v, w (and x)
              // Info(u<<" "<<*iters[0]<<" "<<*iters[1]<<" "<<x<<" "<<y)
            }
          }
        }
      }
      for(auto &w : g.neighOut_iter(u)) parents[w] --;
    }
  }
  Info("Found "<<k<<"-cliques: "<< t << " \t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  return t;
}
ull count_cliques_parallel_edges(const Badjlist &g, const Edgelist &h, ul k, int p) {
  omp_set_num_threads(p);
  Info("Counting cliques for k="<< k<<" with threads "<<p)
  double m = g.e / g.edge_factor;

  ull t = 0, c = 0;
  #pragma omp parallel reduction(+ : t, c)
  {
    vector<ul> parents(g.n, 0);
    #pragma omp for schedule(dynamic, 1)
    for(ul i=0; i<h.e; i++) {
      ul u = h.edges[i].first, v = h.edges[i].second;
      if(g.get_degOut(u) > g.get_degOut(v)) u = h.edges[i].second, v = h.edges[i].first;
      for(auto &w : g.neighOut_iter(u)) parents[w] ++;
      for(auto &w : g.neighOut_iter(v)) parents[w] ++;
      vector<const ul*> iters_end; iters_end.reserve(k-2);
      vector<const ul*> iters; iters.reserve(k-2);
      iters.push_back(g.neighOut_beg(u));
      iters_end.push_back(g.neighOut_end(u));
      while(true) {
        if(iters.back() == iters_end.back()) {
          // Debug("No more neighbours for floor "<< iters.size())
          if(iters.size() == 1) break;
          iters.pop_back();
          iters_end.pop_back();
          for(auto &w : g.neighOut_iter(*iters.back())) parents[w] --;
          ++iters.back();
        }
        else if(iters.size() < k-3) {
          ul v  = *iters.back();
          if(parents[v] < iters.size()+1) { // not all previous nodes are parents
            // Debug("Unsatisfying: "<<parents[v])
            ++iters.back();
            continue;
          }
          for(auto &w : g.neighOut_iter(v)) parents[w] ++;
          iters.push_back(g.neighOut_beg(v));
          iters_end.push_back(g.neighOut_end(v));
        }
        else {
          ul x = *iters.back();
          ++iters.back();
          if(parents[x] < iters.size()+1) continue;
          for(auto &y : g.neighOut_iter(x)) {
            c++;
            if(parents[y] == iters.size()+1) {
              t++; // is neighbour of u, v, w (and x)
              // Info(u<<" "<<*iters[0]<<" "<<*iters[1]<<" "<<x<<" "<<y)
            }
          }
        }
      }
      for(auto &w : g.neighOut_iter(u)) parents[w] --;
      for(auto &w : g.neighOut_iter(v)) parents[w] --;
    }
  }


  Info("Found "<<k<<"-cliques: "<< t << " \t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  return t;
}








//---------------
void set_array(Badjlist &g, vector<ul> &parents, vector<vector<ul>> &deg_level, const ul &u, int l) {
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+deg_level[l-1][u]; ++v) parents[*v] ++;
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+deg_level[l-1][u]; ++v) {
    deg_level[l][*v] = 0;
    for(auto w=g.neigh_beg(*v); w<g.neigh_beg(*v)+deg_level[l-1][*v]; ++w) {
      if(parents[*w] >= l) { // if w is among the neighbours of u, put it at the beginning
        auto x = g.neigh_beg(*v) + deg_level[l][*v]; // neigh with which to swap w
        swap(*w, *x);
        deg_level[l][*v] ++;
      }
    }
  }
}
void reset_array(Badjlist &g, vector<ul> &parents, vector<vector<ul>> &deg_level, const ul &u, int l) {
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+deg_level[l-1][u]; ++v) parents[*v] --;
}

ull count_cliques_parallel_(Badjlist &g, ul k, int p) {
  omp_set_num_threads(p);
  Info("Counting cliques for k="<< k<<" with  "<<p<<"-threads modified")
  double m = g.e / g.edge_factor;

  ull t = 0, c = 0;

  #pragma omp parallel reduction(+ : t, c)
  {
    vector<ul> parents(g.n, 0);
    vector<vector<ul>> deg_level(k-1);
    for(ul l=0; l<k-1; ++l) deg_level[l].reserve(g.n);
    for(ul u=0; u<g.n; ++u) deg_level[0][u] = g.get_degOut(u);

    #pragma omp for schedule(dynamic, 1)
    for(ul u=0; u<g.n; ++u) {
      ul l = 0;
      if(g.get_degOut(u) < k-1) continue;
      set_array(g, parents, deg_level, u, 1);
      vector<const ul*> iters_end; iters_end.reserve(k-1);
      vector<const ul*> iters; iters.reserve(k-1);
      iters.push_back(g.neighOut_beg(u));
      iters_end.push_back(g.neighOut_end(u));
      while(true) {
        l = iters.size();
        if(iters.back() == iters_end.back()) {
          // Debug("No more neighbours for floor "<< l)
          if(l == 1) break;
          iters.pop_back();
          iters_end.pop_back();
          reset_array(g, parents, deg_level, *iters.back(), l);
          ++iters.back();
        }
        else if(l < k-2) {
          ul v  = *iters.back();
          if(deg_level[l][v] < k-1-l) { // not enough out neighbours
            ++iters.back();
            continue;
          }
          set_array(g, parents, deg_level, v, l+1);
          iters.push_back(g.neighOut_beg(v));
          iters_end.push_back(g.neighOut_beg(v)+deg_level[l][v]);
        }
        else {
          ul x = *iters.back();
          ++iters.back();
          for(auto y=g.neigh_beg(x); y<g.neigh_beg(x)+deg_level[l][x]; ++y) {
            c++;
            t++;
          }
        }
      }
      reset_array(g, parents, deg_level, u, 1);
    }
  }
  Info("Found "<<k<<"-cliques: "<< t << " \t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  return t;
}
//---------------

void neigh_subset(Badjlist &g, vector<int> &is_neighOut, vector<ul> &old_deg, vector<ul> &new_deg, ul &u, ul l) {
  // for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+old_deg[u]; ++v) is_neighOut[*v] ++; // set array
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+old_deg[u]; ++v) {
    new_deg[*v] = 0;
    for(auto w=g.neigh_beg(*v); w<g.neigh_beg(*v)+old_deg[*v]; ++w) {
      if(is_neighOut[*w] >= l) { // if w is among the neighbours of u, put it at the beginning
        auto x = g.neigh_beg(*v) + new_deg[*v]; // neigh with which to swap w
        // Info("swap "<<*x<<" "<<*w)
        swap(*w, *x);
        new_deg[*v] ++;
      }
    }
  }
}

void count_rec(Badjlist &g, vector<int> &is_neighOut, vector<vector<ul>> &deg_level, ul &u, ull &t, ull &c, const int k, int l) {
  if(l == k-1) {
    for(auto y=g.neigh_beg(u); y<g.neigh_beg(u)+deg_level[k-2][u]; ++y) {
      c++;
      t++;
    }
    return;
  }
  if(deg_level[l-1][u] < k-l) return;
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+deg_level[l-1][u]; ++v) is_neighOut[*v] ++;
  neigh_subset(g, is_neighOut, deg_level[l-1], deg_level[l], u, l);
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+deg_level[l-1][u]; ++v)
    count_rec(g, is_neighOut, deg_level, *v, t, c, k, l+1);
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+deg_level[l-1][u]; ++v) is_neighOut[*v] --;
}

ull count_cliques_5_(Badjlist &g) {
  Info("Counting 5-cliques with recursion")
  double m = g.e / g.edge_factor;

  ull t = 0, c = 0;
  vector<int> is_neighOut(g.n, 0);
  vector<vector<ul>> deg_level(5-1);
  for(int l=0; l<5-1; ++l) deg_level[l].reserve(g.n);
  for(ul u=0; u<g.n; ++u) deg_level[0][u] = g.get_degOut(u);

  for(ul u=0; u<g.n; ++u) {
    count_rec(g, is_neighOut, deg_level, u, t, c, 5, 1);
  }
  Info("Found 5-cliques: "<< t << " \t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  return t;
}



#define Keep(u, v, deg_level, l, inner) {\
  if(deg_level[l-1][u] < 5-l) continue;\
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+deg_level[l-1][u]; ++v) is_neighOut[*v] ++;\
  neigh_subset(g, is_neighOut, deg_level[l-1], deg_level[l], u, l);\
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+deg_level[l-1][u]; ++v) {\
    inner\
  }\
  for(auto v=g.neigh_beg(u); v<g.neigh_beg(u)+deg_level[l-1][u]; ++v) is_neighOut[*v] --; }

ull count_cliques_5(Badjlist &g) {
  Info("Counting 5-cliques with #macro")
  double m = g.e / g.edge_factor;

  ull t = 0, c = 0;
  vector<int> is_neighOut(g.n, 0);
  vector<vector<ul>> deg_level(5-1);
  for(int l=0; l<5-1; ++l) deg_level[l].reserve(g.n);
  for(ul u=0; u<g.n; ++u) deg_level[0][u] = g.get_degOut(u);

  for(ul u=0; u<g.n; ++u) {
    Keep(u, v, deg_level, 1, {
      Keep(*v, w, deg_level, 2, {
        Keep(*w, x, deg_level, 3, {
          for(auto y=g.neigh_beg(*x); y<g.neigh_beg(*x)+deg_level[3][*x]; ++y) {
            c++;
            t++;
          }
        })
      })
    })
  }
  Info("Found 5-cliques: "<< t << " \t| per edge: " << ((double) t) / m)
  Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
  return t;
}


// ull count_cliques_5_(Badjlist &g) {
//   Info("Counting 5-cliques with modification")
//   double m = g.e / g.edge_factor;
//
//   ull t = 0, c = 0;
//   vector<int> is_neighOut(g.n, 0);
//   vector<ul> deg_lu; deg_lu.reserve(g.n);
//   for(ul u=0; u<g.n; ++u) deg_lu[u] = g.get_degOut(u);
//   vector<ul> deg_lv(g.n, 0);
//   vector<ul> deg_lw(g.n, 0);
//   vector<ul> deg_lx(g.n, 0);
//
//
//   for(ul u=0; u<g.n; ++u) {
//     if(g.get_degOut(u) < 4) continue; // continue if deg u < k-1
//     for(auto v=g.neigh_beg(u); v<g.neigh_end(u); ++v) is_neighOut[*v] ++;
//     neigh_subset(g, is_neighOut, deg_lu, deg_lv, u, 1); // for all v, put its neighbours w at the beginning if they are neigh of u
//     // for(auto v=g.neigh_beg(u); v<g.neigh_end(u); ++v) deg_lv[*v] = deg_lu[*v];
//
//     for(auto v=g.neigh_beg(u); v<g.neigh_end(u); ++v) {
//       if(deg_lv[*v] < 3) continue; // continue if deg v < k-2
//       for(auto w=g.neigh_beg(*v); w<g.neigh_beg(*v)+deg_lv[*v]; ++w) is_neighOut[*w] ++;
//       neigh_subset(g, is_neighOut, deg_lv, deg_lw, *v, 2);
//       // for(auto w=g.neigh_beg(*v); w<g.neigh_beg(*v)+deg_lv[*v]; ++w) deg_lw[*w] = deg_lv[*w];
//
//       for(auto w=g.neigh_beg(*v); w<g.neigh_beg(*v)+deg_lv[*v]; ++w) {
//         if(deg_lw[*w] < 2) continue;
//         for(auto x=g.neigh_beg(*w); x<g.neigh_beg(*w)+deg_lw[*w]; ++x) is_neighOut[*x] ++; // set array
//         neigh_subset(g, is_neighOut, deg_lw, deg_lx, *w, 3);
//
//         for(auto x=g.neigh_beg(*w); x<g.neigh_beg(*w)+deg_lw[*w]; ++x) {
//           for(auto y=g.neigh_beg(*x); y<g.neigh_beg(*x)+deg_lx[*x]; ++y) {
//             c++;
//             t++; // is neighbour of u, v, w (and x)
//             // Info(u<<" "<<v<<" "<<w<<" "<<x<<" "<<y)
//             // if(t > 1000000000) return t;
//           }
//         }
//         for(auto x=g.neigh_beg(*w); x<g.neigh_beg(*w)+deg_lw[*w]; ++x) is_neighOut[*x] --; // reset array
//       }
//       for(auto w=g.neigh_beg(*v); w<g.neigh_beg(*v)+deg_lv[*v]; ++w) is_neighOut[*w] --; // reset array
//     }
//     for(auto v=g.neigh_beg(u); v<g.neigh_end(u); ++v) is_neighOut[*v] --; // reset array
//   }
//   Info("Found 5-cliques: "<< t << " \t| per edge: " << ((double) t) / m)
//   Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
//   return t;
// }


// ull count_cliques_5(const Badjlist &g) {
//   Info("Counting 5-cliques")
//   double m = g.e / g.edge_factor;
//
//   ull t = 0, c = 0;
//   vector<int> is_neighOut(g.n, 0);
//   for(ul u=0; u<g.n; ++u) {
//     for(auto &v : g.neighOut_iter(u)) is_neighOut[v] ++; // set array
//     for(auto &v : g.neighOut_iter(u)) {
//       for(auto &w : g.neighOut_iter(v)) is_neighOut[w] ++; // set array
//       for(auto &w : g.neighOut_iter(v)) {
//         if(is_neighOut[w] < 2) continue;
//         for(auto &x : g.neighOut_iter(w)) is_neighOut[x] ++; // set array
//         for(auto &x : g.neighOut_iter(w)) {
//           if(is_neighOut[x] < 3) continue;
//           for(auto &y : g.neighOut_iter(x)) {
//             c++;
//             if(is_neighOut[y] == 3) t++; // Info(u<<" "<<v<<" "<<w<<" "<<x<<" "<<y)
//           }
//         }
//         for(auto &x : g.neighOut_iter(w)) is_neighOut[x] --; // reset array
//       }
//       for(auto &w : g.neighOut_iter(v)) is_neighOut[w] --; // reset array
//     }
//     for(auto &v : g.neighOut_iter(u)) is_neighOut[v] --; // reset array
//   }
//   Info("Found 5-cliques: "<< t << " \t| per edge: " << ((double) t) / m)
//   Info("Complexity: "<< c << " \t| per edge: " << ((double) c) / m)
//   return t;
// }
