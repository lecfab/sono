// See header for general documentation

#include <functional>
#include <unordered_set>
#include <set>
#include "order_triangles.h"
#include "../utils/adjlist.h"
#include "../utils/continuousrank.h"

using namespace std;

// mistakes: pi should have double type in neighbour_positions
// mistakes: end value should be pi+1 and not pointer+1
// improved: watchlist of updated nodes
// improved: linear computing of best penalties
// improved: linked list instead of array where everyone had to be moved one by one
// improved: when u goes after v, put u close to v instead of at the end of the chain (normalisation and cache purpose)
// improved: ContinuousRank is now a separate structure with hidden tools
// todo: watch heap instead of bool-array, but it seems slower in practice
// todo: altitude to be evened out among several predecessors
vector<ul> place_neighbour_dpp(const Adjlist &g)    { return place_neighbour(g, sum_square_dpp, best_position_dpp); }
vector<ul> place_neighbour_dpm(const Adjlist &g)    { return place_neighbour(g, sum_square_dpm, best_position_dpm); }
vector<ul> place_neighbour_dpm2(const Adjlist &g)    { return place_neighbour(g, sum_dpm2, best_position_dpm2); }
vector<ul> place_neighbour_mindpp(const Adjlist &g) { return place_neighbour(g, sum_mindpp,  best_position_mindpp); }
vector<ul> place_neighbour(const Adjlist &g, decltype(sum_square_dpm) &scoring, decltype(best_position_dpm) &best_position) {
  ul m = g.e / g.edge_factor;
  // ----- preparation -----
  ContinuousRank cRank(g.n);
  // unordered_set<ul> next_nodes; next_nodes.reserve(g.n);

  vector<ul> dp = compute_dp(g, cRank);
  vector<bool> watch(g.n, true);

  // ----- loop until stability is reached -----
  ul score = scoring(g, dp), previous_score = score * 2;
  cout << "Initial score: " << ((double) score) / m << endl;
  double precision = 1*1e-3; // SET PRECISION
  while(score < previous_score) {
    if(score >= (1.-precision) * previous_score) {
      Info("Truncated with precision "<< precision)
      break;
    }
    previous_score = score;

    // select nodes that have been updated
    vector<ul> nodes; nodes.reserve(g.n);
    for(ul u=0; u<g.n; ++u)
      if(watch[u]) {
        watch[u] = false;
        nodes.push_back(u);
      }

    // ----- consider each node in random order -----
    random_shuffle(nodes.begin(), nodes.end());
    for(auto &u : nodes) {
      vector<ul> neighbours; neighbours.reserve(g.get_degree(u));
      for(auto &v : g.neigh_iter(u)) neighbours.push_back(v);
      sort(neighbours.begin(), neighbours.end(), [&cRank](ul v, ul w) { return cRank[v] < cRank[w]; });
      // ----- find the best position in its neighbourhood -----
      ul i0 = g.get_degree(u) - dp[u];
      auto pos_penalty = best_position(g, u, neighbours, dp, i0);
      if(pos_penalty.second == 0) continue;

      // ----- put it in the best position: update score, out-degrees and rank -----
      if(pos_penalty.second > 0) // code-word to say that dp has already been updated
        score -= pos_penalty.second;
      else {
        score += pos_penalty.second;
        dp[u] = g.get_degree(u) - pos_penalty.first;
        for(ul j=pos_penalty.first; j<i0; ++j) dp[ neighbours[j] ] --;
        for(ul j=i0; j<pos_penalty.first; ++j) dp[ neighbours[j] ] ++;
      }

      if(pos_penalty.first == g.get_degree(u)) cRank.move_after(u, neighbours.back());
      else cRank.move_before(u, neighbours[pos_penalty.first]);

      // for(auto &v : neighbours) { Info("\t"<< cRank[v])}

      for(auto &v : neighbours) watch[v] = true;
    }

    cout << "Current score: " << ((double) score) / m << "\t ("<<100*((double) nodes.size()) / g.n<<"% nodes) \tNormalisations: "<<cRank.normalisations <<endl;
  }

  // just to check that score is accurate
  // dp = compute_dp(g, cRank);
  // score = scoring(g, dp);
  cout << "Final score:   " << ((double) score) / m << endl;

  return cRank.to_rank();
}

vector<ul> compute_dp(const Adjlist &g, const ContinuousRank &cRank) {
  vector<ul> dp(g.n, 0);
  for(ul u=0; u<g.n; ++u)
    for(auto &v : g.neigh_iter(u))
      if(cRank[v] > cRank[u]) dp[u] ++;
  return dp;
}



ul sum_square_dpp(const Adjlist &g, const vector<ul> &vec) {
  ul s = 0;
  for(auto &a : vec)
    s += pow(a, 2);
  return s;
}
ul sum_square_dpm(const Adjlist &g, const vector<ul> &dp) {
  ul s = 0;
  for(ul u=0; u<g.n; ++u)
    s += dp[u] * (g.get_degree(u) - dp[u]);
  return s;
}
ul sum_dpm2(const Adjlist &g, const vector<ul> &dp) { return sum_dpm_ab(g, dp, 2, 2); }
ul sum_dpm_ab(const Adjlist &g, const vector<ul> &dp, double a, double b) {
  ul s = 0;
  for(ul u=0; u<g.n; ++u)
    s += pow(dp[u], a) * pow(g.get_degree(u) - dp[u], b);
  return s;
}
ul sum_mindpp(const Adjlist &g, const vector<ul> &dp) { // , const vector<double> &altitude
  ul s = 0;
  for(ul u=0; u<g.n; ++u)
    for(auto &v : g.neigh_iter(u))
      if(v > u) // if(altitude[v] > altitude[u]) would be more correct but we just want to select half of the edges
        s += min(dp[u], dp[v]);
  return s;
}




pair<ul, long long int> best_position_mindpp(const Adjlist &g, const ul &u, const vector<ul> &neighbours, vector<ul> &dp, ul i0) {
  pair<ul, long long int> pos_penalty = { 0, 0 };
  for(ul i=i0+1; i<=g.get_degree(u); ++i) {
    ul v = neighbours[i-1];
    // ----- compute the variation if u goes after v -----
    long long int delta = (dp[u] > dp[v] + 1);
    for(auto &w : g.neigh_iter(u)) if(dp[u] <= dp[w]) delta --;
    for(auto &w : g.neigh_iter(v)) if(w != u and dp[v] < dp[w]) delta ++;

    if(delta >= 0) break;
    // ----- update score and out-degrees -----
    pos_penalty = { i, pos_penalty.second + delta };
    dp[u] --;
    dp[v] ++;
  }
  return { pos_penalty.first, -pos_penalty.second }; // second is inverted: code-word to say that dp has already been updated
}

pair<ul, long long int> best_position_dpp(const Adjlist &g, const ul &u, const vector<ul> &neighbours, vector<ul> &dp, ul i0) {
  pair<ul, long long int> pos_penalty = { 0, 0 };
  long long int neighbour_penalty = 0;
  for(ul i=i0; i>0; --i) {
    ul j = i-1;
    neighbour_penalty += 1 - 2*dp[neighbours[j]];
    long long int penalty = neighbour_penalty + j*j - i0*i0 + 2*(i0-j)*g.get_degree(u);
    if(penalty < pos_penalty.second) pos_penalty = { j, penalty };
  }
  neighbour_penalty = 0;
  for(ul i=i0; i<g.get_degree(u); ++i) {
    ul j = i+1;
    neighbour_penalty += 1 + 2*dp[neighbours[i]];
    long long int penalty = neighbour_penalty + j*j - i0*i0 + 2*(i0-j)*g.get_degree(u);
    if(penalty < pos_penalty.second) pos_penalty = { j, penalty };
  }
  return pos_penalty;
}

pair<ul, long long int> best_position_dpm(const Adjlist &g, const ul &u, const vector<ul> &neighbours, vector<ul> &dp, ul i0) {
  pair<ul, long long int> pos_penalty = { 0, 0 };
  long long int neighbour_penalty = 0;
  for(ul i=i0; i>0; --i) {
    ul j = i-1;
    neighbour_penalty += 2*dp[neighbours[j]] - g.get_degree(neighbours[j]) - 1;
    long long int penalty = neighbour_penalty + i0*i0 - j*j + (j-i0)*g.get_degree(u);
    if(penalty < pos_penalty.second) pos_penalty = { j, penalty };
  }
  neighbour_penalty = 0;
  for(ul i=i0; i<g.get_degree(u); ++i) {
    ul j = i+1;
    neighbour_penalty += g.get_degree(neighbours[i]) - 2*dp[neighbours[i]] - 1;
    long long int penalty = neighbour_penalty + i0*i0 - j*j + (j-i0)*g.get_degree(u);
    if(penalty < pos_penalty.second) pos_penalty = { j, penalty };
  }

  // auto pp_ab = best_position_dpm_ab(g, u, neighbours, dp, i0, 1, 1);
  // if(pp_ab.second != pos_penalty.second) {
  //   Alert(pp_ab.second <<"!="<< pos_penalty.second)
  //   exit(0);
  // }
  return pos_penalty;
}


pair<ul, long long int> best_position_dpm2(const Adjlist &g, const ul &u, const vector<ul> &neighbours, vector<ul> &dp, ul i0) {
  return best_position_dpm_ab(g, u, neighbours, dp, i0, 2, 2);
}

pair<ul, long long int> best_position_dpm_ab(const Adjlist &g, const ul &u, const vector<ul> &neighbours, vector<ul> &dp, ul i0, double a, double b) {
  pair<ul, long long int> pos_penalty = { 0, 0 };
  long long int neighbour_penalty = 0;
  for(ul i=i0; i>0; --i) {
    ul j = i-1;
    ul v = neighbours[j];
    neighbour_penalty += pow(dp[v]-1, a)*pow(g.get_degree(v)-dp[v]+1, b) - pow(dp[v], a)*pow(g.get_degree(v)-dp[v], b);
    long long int penalty = neighbour_penalty + pow(g.get_degree(u)-j, a)*pow(j, b) - pow(dp[u], a)*pow(g.get_degree(u)-dp[u], b);
    if(penalty < pos_penalty.second) pos_penalty = { j, penalty };
  }
  neighbour_penalty = 0;
  for(ul i=i0; i<g.get_degree(u); ++i) {
    ul j = i+1;
    ul v = neighbours[i]; // index i is valid (rank in the neighbourhood), not j (rank if u is also counted)
    neighbour_penalty += pow(dp[v]+1, a)*pow(g.get_degree(v)-dp[v]-1, b) - pow(dp[v], a)*pow(g.get_degree(v)-dp[v], b);
    long long int penalty = neighbour_penalty + pow(g.get_degree(u)-j, a)*pow(j, b) - pow(dp[u], a)*pow(g.get_degree(u)-dp[u], b);
    if(penalty < pos_penalty.second) pos_penalty = { j, penalty };
  }
  return pos_penalty;
}





// pair<ul, long long int> best_position_mindpp(const Adjlist &g, const ul &u, const vector<ul> &neighbours, const vector<ul> &dp, ul i0) {
//   pair<ul, long long int> pos_penalty = { 0, 0 };
//   ul new_dp_u = 0;
//   vector<ul> new_dp; new_dp.reserve(neighbours.size());
//   for(ul j=0; j<i0; ++j) new_dp.push_back(dp[neighbours[j]] - 1); // predecessors loose a successor if u is in position 0
//   for(ul j=i0; j<g.get_degree(u); ++j) new_dp.push_back(dp[neighbours[j]]); // successors are not impacted if u is in position 0
//
//   for(ul i=0; i<=g.get_degree(u); ++i) { // check each new position for u
//     new_dp_u ++; // one more predecessor for u
//     if(i > 0) new_dp[i-1] ++; // one more successor for this predecessor
//
//     long long int penalty = 0;
//     for(auto &v : g.neigh_iter(u)) penalty += min(new_dp_u, new_dp[v]) - min(dp[u], dp[v]);
//     ...
//
//
//     long long int penalty = i0*i0 - i*i + (i-i0)*g.get_degree(u);
//     for(ul j=i; j<i0; ++j) // move neighbours to the right if i<i0
//       penalty += 2*dp[neighbours[j]] - g.get_degree(neighbours[j]) - 1;
//     for(ul j=i0; j<i; ++j) // move neighbours to the left if i>i0
//       penalty += g.get_degree(neighbours[j]) - 2*dp[neighbours[j]] - 1;
//     if(penalty < best_penalty) {
//       best_pos = i; best_penalty = penalty;
//     }
//   }
// }


// vector<ul> optimise_edges(const Adjlist &g) {
//   ul m = g.e / g.edge_factor;
//   // ----- preparation -----
//   vector<double> altitude; altitude.reserve(g.n);
//   list<ul> chain;
//   vector<list<ul>::iterator> pointer; pointer.reserve(g.n);
//   for(ul u=0; u<g.n; ++u) {
//     altitude.push_back(u-((double) g.n) / 2);
//     pointer.push_back(chain.insert(chain.end(), u));
//   }
//
//   vector<ul> dp = compute_dp(g, altitude);
//   vector<bool> watch(g.n, true);
//
//   // ----- loop until stability is reached -----
//   ul score = sum_mindpp(g, dp), previous_score = score * 2;
//   cout << "Initial score: " << ((double) score) / m << endl;
//   double precision = 0*1e-3; // SET PRECISION
//   while(score < previous_score) {
//     if(score >= (1.-precision) * previous_score) {
//       Info("Truncated with precision "<< precision)
//       break;
//     }
//     previous_score = score;
//
//     // select nodes that have been updated
//     vector<ul> nodes; nodes.reserve(g.n);
//     for(ul u=0; u<g.n; ++u)
//       if(watch[u]) {
//         watch[u] = false;
//         nodes.push_back(u);
//       }
//
//     // ----- consider each node in random order -----
//     random_shuffle(nodes.begin(), nodes.end());
//     for(auto &u : nodes) {
//       vector<ul> neighbours; neighbours.reserve(g.get_degree(u));
//       for(auto &v : g.neigh_iter(u)) if(altitude[v] > altitude[u]) neighbours.push_back(v);
//       sort(neighbours.begin(), neighbours.end(), [&altitude](ul v, ul w) { return altitude[v] < altitude[w]; });
//
//       // ----- select the closest successor according to current ordering -----
//       ul v_final = u;
//       for(auto &v : neighbours) {
//         // ----- compute the variation if u goes after v -----
//         long long int delta = (dp[u] > dp[v] + 1);
//         for(auto &w : g.neigh_iter(u)) if(dp[u] <= dp[w]) delta --;
//         for(auto &w : g.neigh_iter(v)) if(w != u and dp[v] < dp[w]) delta ++;
//
//         if(delta >= 0) break;
//         // ----- update score and out-degrees -----
//         score += delta;
//         dp[u] --;
//         dp[v] ++;
//
//         v_final = v;
//       }
//       if(v_final == u) continue;
//
//       // ----- put it in the new position: update altitude and watchlist -----
//       pointer[u] = insert_into_chain(chain, pointer[u], next(pointer[v_final]), altitude);
//       watch[u] = true;
//       for(auto &w : g.neigh_iter(u)) // new predecessors of u should be watched
//         if(altitude[w] < altitude[u]) watch[w] = true;
//     }
//
//     cout << "Current score: " << ((double) score) / m << "\t ("<<100*((double) nodes.size()) / g.n<<"% nodes)" <<endl;
//   }
//
//   // just to check that score is accurate
//   dp = compute_dp(g, altitude);
//   score = sum_mindpp(g, dp);
//   cout << "Final score:   " << ((double) score) / m << endl;
//
//   return rank_from_chain(chain);
// }
