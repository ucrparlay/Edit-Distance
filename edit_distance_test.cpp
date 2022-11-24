#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include "edit_distance_dp.h"
#include "edit_distance_parallel.h"
#include "edit_distance_sequential.h"
#include "parlay/internal/get_time.h"

int main(int argc, char** argv) {
  srand(time(0));
  EditDistanceDP ed;
  EditDistanceParallel ed_par;
  parlay::internal::timer t1("t1", false), t2("t2", false);
  for (int tt = 0; tt < 1; tt++) {
    // int n = rand() % 10000 + 1;
    // int m = rand() % 10000 + 1;
    int n = 100000, m = 10000;
    std::string s1, s2;
    for (int i = 0; i < n; i++) {
      s1 += 'a' + rand() % 3;
    }
    for (int i = 0; i < m; i++) {
      s2 += 'a' + rand() % 3;
    }
    // s1 = "cbcbcbb", s2 = "c";
    t1.start();
    int d1 = ed.Solve(s1, s2);
    t1.stop();
    t2.start();
    int d2 = ed_par.Solve(s1, s2);
    t2.stop();
    if (d1 != d2) {
      std::cout << d1 << ' ' << d2 << std::endl;
      std::cout << s1 << '\n' << s2 << '\n' << n << '\n' << m << std::endl;
      assert(d1 == d2);
    }
  }
  std::cout << "Edit Distance Test Pass" << std::endl;
  t1.total();
  t2.total();
}