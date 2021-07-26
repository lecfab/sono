#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <fstream>      // ofstream: write, ifstream: read, fstream: read and write from/to files
#include <algorithm>    // std::random_shuffle
#include <ctime>        // std::time
// #include <time.h>       // time_t, struct tm, difftime, time, mktime
#include <cstdlib>      // std::rand, std::srand
#include <string>
#include <chrono>
#include <cmath>
#include <vector>
#include <utility>
#include <list>

#ifndef DEBUG
#define Debug(msg) {}
#else
#define Debug(msg) { std::cout << "\tDebug: " << msg << endl; fflush(stdout); }
#endif
#define Alert(msg) { std::cout << "\tAlert: " << msg << endl; fflush(stdout); }
#define Info(msg) { std::cout << "\tInfo: " << msg << endl; fflush(stdout); }

#ifndef NLINKS
#define NLINKS 100000000
#define NNODES 10000000
#endif

#define TimeBegin() auto time_begin = std::chrono::steady_clock::now(); \
                    auto time_previous = time_begin, timeRec_previous = time_begin, time_now = time_begin; \
                    std::time_t t0, t1, t2; t1 = std::time(nullptr); t0 = t1;
                    // auto clock_begin = std::clock();
                    // auto clock_previous = clock_begin, clockRec_previous = clock_begin, clock_now = clock_begin;
                    // auto time_diff = std::chrono::duration_cast<std::chrono::milliseconds>(time_now - time_previous).count();
                    // auto clock_diff = clock_now - clock_begin;
#define TimeStep(m) std::cout << m; \
                    time_now = std::chrono::steady_clock::now(); \
                    t2 = time(NULL); \
                    print_chrono_ctime(time_previous, time_now, t1, t2); \
                    time_previous = time_now; \
                    t1 = t2;
                    // clock_now = std::clock();
                    // clock_previous = clock_now;
                    // print_chrono_clock(time_previous, time_now, clock_previous, clock_now);
#define TimeTotal()	std::cout << "Total"; \
                    time_now = std::chrono::steady_clock::now(); \
                    t2 = time(NULL); \
                    print_chrono_ctime(time_begin, time_now, t0, t2);
                    // clock_now = std::clock();
                    // print_chrono_clock(time_begin, time_now, clock_begin, clock_now);

#define TimeRecStep(m, out) time_now = std::chrono::steady_clock::now(); \
                            t2 = time(NULL); \
                            out << m << "\t"; \
                            out << std::chrono::duration_cast<std::chrono::milliseconds>(time_now - timeRec_previous).count() << "\t"; \
                            out << t2 - t1 << endl; \
                            timeRec_previous = time_now; \
                            t1 = t2;
                            // clock_now = std::clock();
                            // clockRec_previous = clock_now;
                            // out << clock_now - clockRec_previous << endl;
#define TimeRecTotal(m, out)  time_now = std::chrono::steady_clock::now(); \
                              t2 = time(NULL); \
                              out << m << "\t"; \
                              out << std::chrono::duration_cast<std::chrono::milliseconds>(time_now - time_begin).count() << "\t"; \
                              out << t2 - t0 << endl;
                              // clock_now = std::clock();
                              // out << clock_now - clock_begin << endl;

// #define TimeRecStep(m) time_now = chrono::steady_clock::now(); std::cout << chrono::duration_cast<chrono::milliseconds>(time_now - time_previous).count() << "\t"; time_previous = time_now; Debug(m)
// #define TimeRecTotal() time_now = chrono::steady_clock::now(); std::cout << chrono::duration_cast<chrono::milliseconds>(time_now - time_begin).count() << endl;

// #define TimeRecBegin() vector<chrono::steady_clock::time_point> time_rec; time_rec.push_back(chrono::steady_clock::now());
// #define TimeRecStep(m) time_rec.push_back(chrono::steady_clock::now());
// #define TimeRecTotal() { auto t0 = time_rec[0], t1=t0; for(auto &t : time_rec) { std::cout << chrono::duration_cast<chrono::milliseconds>(t - t1).count() << ","; t1=t; } std::cout << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << endl; }

void print_chrono_ctime(
  std::chrono::steady_clock::time_point begin,
  std::chrono::steady_clock::time_point end,
  std::time_t t_begin,
  std::time_t t_end);
void print_chrono_clock(
  std::chrono::steady_clock::time_point begin,
  std::chrono::steady_clock::time_point end,
  std::clock_t cl_begin,
  std::clock_t cl_end);
void print_chrono(std::chrono::steady_clock::time_point begin, std::chrono::steady_clock::time_point end);

// https://www.cprogramming.com/c++11/c++11-auto-decltype-return-value-after-function.html
// template <typename Builder>
// auto makeAndProcessObject (const Builder& builder) -> decltype( builder.makeObject() );
// template <typename T>
// auto Timer (T &function, int argc, char** argv) -> decltype(function(argc, argv));


typedef unsigned long ul;
typedef unsigned long long ull;
//compute the maximum of three unsigned long
inline ul max3(ul a, ul b, ul c) { return (a >= b and a >= c) ? a : ((b >= c) ? b : c); }

struct Keyvalue {
	ul key;
	ul val;
	Keyvalue(ul key, ul val) : key(key), val(val) {}
};

std::vector<ul> order_identity(const ul &n);
std::vector<ul> rank_from_order(const std::vector<ul> &order, const ul &n);
std::vector<ul> rank_from_order(const std::vector<ul> &order);
std::pair<ul,ul> random_pair(const ul &n);

#endif
