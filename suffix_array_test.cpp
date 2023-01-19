#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "parlay/internal/get_time.h"
#include "parlay/parallel.h"
#include "suffix_array_parallel.h"
#include "suffix_array_sequential.h"
#include "hash_lcp_parallel.h"

using namespace std;

int main() {
  int n = 1000, m = 1000;
  parlay::sequence<int> a(n), b(m);
  for (int i = 0; i < n; i++) {
    a[i] = rand() % 2;
  }
  for (int i = 0; i < m; i++) {
    b[i] = rand() % 2;
  }
  // parlay::internal::timer timer1;
  // auto [rank1, sa1, lcp1] = suffix_array_large_alphabet(s);
  // auto t1 = timer1.stop();
  // cout << "Parallel time: " << to_string(t1) << endl;
  // parlay::internal::timer timer2;
  // SuffixArraySequential sa_seq;
  // sa_seq.Build(s);
  // auto t2 = timer2.stop();
  // auto rank2 = sa_seq.rank;
  // auto sa2 = sa_seq.sa;
  // auto lcp2 = sa_seq.height;
  // assert(sa1.size() == sa2.size());
  // for (int i = 0; i < n; i++) {
  //   assert(rank1[i] == rank2[i]);
  //   assert(sa1[i] == sa2[i]);
  //   assert(lcp1[i] == lcp2[i]);
  // }
  // cout << "Sequential time: " << to_string(t2) << endl;
  // cout << "Pass: suffix_array_test" << endl;

  vector<int> logN1;
  vector<int> powerN1;
  vector<vector<int>> table_s1;
  vector<vector<int>> table_s2;
  build_hash_table(a, b, table_s1, table_s2, powerN1, logN1);
  std::cout << "hash build done" << std::endl;
  query_lcp(a, b, table_s1, table_s2, logN1, 0, 0);


  auto c = parlay::sequence<uint32_t>(n + m + 1);
  parlay::parallel_for(0, n, [&](int i) { c[i] = a[i]; });
  parlay::parallel_for(0, m, [&](int i) { c[i + n + 1] = b[i]; });
  c[n] = std::numeric_limits<uint32_t>::max();
  auto rank = parlay::sequence<unsigned int>();
  auto sa = parlay::sequence<unsigned int>();
  auto lcp = parlay::sequence<unsigned int>();
  std::tie(rank, sa, lcp) = suffix_array_large_alphabet(c);
  auto rmq = range_min(lcp);
  auto GetLcp = [&](int i, int j) -> int {
    // std::cout << "GetLcp " << i << ' ' << j << '\n';
    if (i == n || j == m) return 0;
    if (a[i] != b[j]) return 0;
    int l = rank[i], r = rank[j + n + 1];
    if (l > r) std::swap(l, r);
    int id = rmq.query(l + 1, r);
    return lcp[id];
  };

  int q = 10000000;
  vector<pair<int, int>> queries(q);
  for (int i = 0; i < q; i++) {
    int x = rand() % n;
    int y = rand() % m;
    queries[i] = {x, y};
  }

  for (auto [i, j] : queries) {
    int x = query_lcp(a, b, table_s1, table_s2, logN1, i, j);
    int y = GetLcp(i, j);
    assert(x == y);
  }

  parlay::internal::timer t2;
  for (auto [i, j] : queries) {
    GetLcp(i, j);
  }
  cout << t2.stop() << std::endl;

  parlay::internal::timer t1;
  for (auto [i, j] : queries) {
    query_lcp(a, b, table_s1, table_s2, logN1, i, j);
  }
  cout << t1.stop() << std::endl;

  return 0;
}
