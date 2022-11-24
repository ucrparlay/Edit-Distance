#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "parlay/internal/get_time.h"
#include "parlay/parallel.h"
#include "suffix_array_parallel.h"
#include "suffix_array_sequential.h"

using namespace std;

int main() {
  int n = 1000000;
  string s;
  s.resize(n);
  for (int i = 0; i < n; i++) {
    s[i] = 'a' + rand() % 2;
  }
  parlay::internal::timer timer1;
  auto [rank1, sa1, lcp1] = suffix_array(s);
  auto t1 = timer1.stop();
  cout << "Parallel time: " << to_string(t1) << endl;
  parlay::internal::timer timer2;
  SuffixArraySequential sa_seq;
  sa_seq.Build(s);
  auto t2 = timer2.stop();
  auto rank2 = sa_seq.rank;
  auto sa2 = sa_seq.sa;
  auto lcp2 = sa_seq.height;
  assert(sa1.size() == sa2.size());
  for (int i = 0; i < n; i++) {
    assert(rank1[i] == rank2[i]);
    assert(sa1[i] == sa2[i]);
    assert(lcp1[i] == lcp2[i]);
  }
  cout << "Sequential time: " << to_string(t2) << endl;
  cout << "Pass: suffix_array_test" << endl;
  return 0;
}
