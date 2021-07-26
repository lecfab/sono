// http://www.cplusplus.com/forum/articles/10627/ to see how to include headers

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "utils/CLI11.h" // options parser
#include "utils/tools.h"
#include "utils/inout.h"
#include "utils/edgelist.h"
#include "utils/adjlist.h"
#include "order/order_rand.h"
#include "order/order_deg.h"
#include "algo/algo_nq.h"
#include "algo/algo_bfs.h"
#include "algo/algo_dfs.h"
#include "algo/algo_bellman.h"
#include "algo/algo_diameter.h"
#include "algo/algo_kcore.h"
#include "algo/algo_tarjan.h"
#include "algo/algo_pagerank.h"
#include "algo/algo_dominatingset.h"
#include "algo/algo_minuncut.h"
#include "algo/algo_triangles.h"
#include "order/order_triangles.h"

using namespace std;


int main(int argc, char** argv){
	TimeBegin()

  CLI::App app{"Execute various algorithms on a given graph."};

  string filename, algo_name="no", output_file="";
  bool directed, algo_all=false; int attempts=1; int number=0;
  app.add_option("file", filename, "Text file: list of `a b` edges with nodes IDs ranging from 0 to N-1")->required();
  app.add_flag("-d,!-u,--directed,!--undirected", directed, "Specify if the graph is directed or undirected; multiple edges are not accepted");
  app.add_option("-a,--algo", algo_name, "Algorithm to execute on the graph")
    ->check(CLI::IsMember({
      "nq",
      "bfs","dfs","tarjan",
      "bellman","diameter",
      "kcore","ds","dominatingset",
      "pr",
      "uncut","triangles","trianglesdpp","trianglesdpm","cliques",
      "test"}, CLI::ignore_case));
  // app.add_flag("-A,--all", algo_all, "Execute all algorithms consecutively and save time measurements")->capture_default_str();
  app.add_option("-o,--output", output_file, "File in which to output time measurements")->capture_default_str();
  app.add_option("-l,--loops", attempts, "Number of attempts for random algorithms")->capture_default_str();
  app.add_option("-n", number, "Numeric parameter used in different orders")->capture_default_str();

  CLI11_PARSE(app, argc, argv);
  ofstream output(output_file);
  // streambuf *buf; ofstream of;
  // if(output_file != "") {
  //   Info("Measurements will be saved in " << output_file)
  //   of.open(output_file);
  //   buf = of.rdbuf();
  // }
  // else {
  //   Info("Measurements will be printed hereafter")
  //   buf = cout.rdbuf();
  // }
  // ostream output(buf);

	// --------------------------------------------------
	// -- Convert edgelist file into Adjlist structure --
	// --------------------------------------------------
	Info("Reading edgelist from file " << filename)

	ifstream file(filename);
	Edgelist h(file); //*h = c_readedgelist(argv[1]);

	Info("Number of nodes: " << h.n)
	Info("Number of edges: " << h.e)

	TimeStep("Read")
  TimeRecStep("read", output)

  srand (time(NULL));

  if(algo_name == "triangles" or algo_name == "cliques") { // uses Badjlist
    Info("Converting to bidirected adjacency list")
    // vector<ul> rank = order_deg(h, /*DESC=*/true);
    // vector<ul> rank = order_identity(h.n);
    // random_shuffle( rank.begin(), rank.end() );
    // h.apply_rank(rank);
    h.to_dag();

    Badjlist g(h);
    TimeStep("Adjlist")
    TimeRecStep("adjlist", output)

    Info("Neighbourhoods are " << (g.are_neighbours_sorted() ? "":"not ") << "sorted")
    // algo_kcore(g); TimeStep("kcore")
    triangle_complexities(g); TimeStep("t0")

    if(algo_name == "triangles") {
      // count_triangles_dpp(g); TimeStep("t1") TimeRecStep("t1", output)
      // count_triangles_dmm(g); TimeStep("t2") TimeRecStep("t2", output)
      count_triangles_dpm(g); TimeStep("t3") TimeRecStep("t3", output)
      count_triangles_dichotomy(g); TimeStep("t4") TimeRecStep("t4", output)
      count_cliques(g, 3); TimeStep("k3") TimeRecStep("k3", output)
      // count_triangles_dpm_indep_vectint(g); TimeStep("t6") TimeRecStep("t6", output)
      // count_triangles_dpm_indep(g); TimeStep("t5") TimeRecStep("t5", output)
      // count_triangles_hash(g); TimeStep("t5")
    }
    else {
      if(number == 0) number = 5;
      if(number == 5) { count_cliques_5(g); TimeStep("k5") TimeRecStep("k5", output) }
      count_cliques(g, number); TimeStep("ki") TimeRecStep("ki", output)
      count_cliques_parallel(g, number, 8); TimeStep("pn8") TimeRecStep("pn8", output)
      count_cliques_parallel_edges(g, h, number, 8); TimeStep("pe8") TimeRecStep("pe8", output)
      count_cliques_parallel(g, number, 16); TimeStep("pn16") TimeRecStep("pn16", output)
      count_cliques_parallel_edges(g, h, number, 16); TimeStep("pe16") TimeRecStep("pe16", output)
    }
  }

  else if(algo_name == "uncut" and !directed) {
    algo_minuncut(h);
    TimeStep("Uncut")
    TimeRecStep("uncut", output)
  }

  else if(algo_name == "pr") { // uses Edgelist and not Adjlist
    algo_pagerank(h, attempts*attempts); // attempts
    TimeStep("PR")
    TimeRecStep("pr", output)
  }

  else { // requires Adjlist
    Info("Converting to " << (directed ? "":"un") << "directed adjacency list")
    Adjlist* g;
    if(directed) g = new Dadjlist(h);
    else g = new Uadjlist(h);

  	TimeStep("Adjlist")
    TimeRecStep("adjlist", output)

    if(algo_all or algo_name != "no") {
  		// --------------------------------------------------
  		// ----- Execute an algorithm over the Adjlist ------
  		// --------------------------------------------------
      TimeTotal()
      TimeRecTotal("prepare", output)


      if(algo_name == "trianglesdpp" and !directed) {
        place_neighbour_dpp(*g);
        TimeStep("trianglesdpp")
        TimeRecStep("trianglesdpp", output)
      }

      if(algo_name == "trianglesdpm" and !directed) {
        place_neighbour_dpm(*g);
        TimeStep("trianglesdpm")
        TimeRecStep("trianglesdpm", output)
      }

  		if(algo_all or algo_name == "nq") {
        algo_nq(*g, g->n);
        TimeStep("NQ")
        TimeRecStep("nq", output)
      }
  		if(algo_all or algo_name == "bfs") {
        for (int i = 0; i < attempts; i++) algo_bfs(*g, rand() % g->n);
        TimeStep("BFS")
        TimeRecStep("bfs", output)
      }
  		if(algo_all or algo_name == "dfs") {
        for (int i = 0; i < attempts; i++) algo_dfs(*g, rand() % g->n);
        TimeStep("DFS")
        TimeRecStep("dfs", output)
      }
      if(algo_all or algo_name == "tarjan") {
        for (int i = 0; i < attempts; i++) algo_tarjan(*g, rand() % g->n);
        TimeStep("SCC")
        TimeRecStep("tarjan", output)
      }
      if(algo_all or algo_name == "bellman") {
        for (int i = 0; i < attempts; i++) algo_bellman(*g, rand() % g->n);
        TimeStep("SP")
        TimeRecStep("bellman", output)
      }

      if(algo_all or algo_name == "ds" or algo_name == "dominatingset") {
        algo_dominatingset(*g);
        TimeStep("DS")
        TimeRecStep("dominatingset", output)
      }
      if(algo_all or algo_name == "kcore") {
        algo_kcore(*g);
        TimeStep("Kcore")
        TimeRecStep("kcore", output)
      }
      if(algo_all or algo_name == "diameter") {
        algo_diameter(*g, attempts * attempts);
        TimeStep("Diam")
        TimeRecStep("diameter", output)
      }
      // -----------------------------------------
  	}

    // delete g; // delete h; // see unique_ptr instead
  }

	TimeTotal()
  TimeRecTotal("total", output);
	return 0;
}
