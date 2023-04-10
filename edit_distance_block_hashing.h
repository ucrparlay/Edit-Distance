#include "hash_block_parallel.h"

// Returns an 32-bit integer of edit distance for two sequences A and B
// using bfs and hashing table query for lcp, in parallel
template <typename Seq>
int EditDistanceBlockHashParallel(const Seq &a, const Seq &b,
                                  double *building_tm) {
  parlay::internal::timer tt;
  // build sparse table

  /**
   * rocord of building time
   */
  Timer tmr;

  int n = a.size();
  int m = b.size();
  if (n == 0) return m;
  if (m == 0) return n;
  parlay::sequence<parlay::sequence<int>> table_A;
  parlay::sequence<parlay::sequence<int>> table_B;
  parlay::sequence<std::pair<int, int>> S_A;
  parlay::sequence<std::pair<int, int>> S_B;
  construct_table(a, b, table_A, table_B, S_A, S_B, std::min(n, m));
  auto Diag = [&](int i, int j) { return i - j + m; };
  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);

  *building_tm = tmr.elapsed();
  tt.next("building");
  // std::cout << building_tm << ", ";

  max_row[Diag(0, 0)] = block_query_lcp(0, 0, a, b, table_A, table_B, S_A, S_B);
  // bfs for path
  int k = 0;
  int round = 0;
  for (;;) {
    // printf("round: %d\n", round++);
    // tt.next("Query");
    if (max_row[Diag(n, m)] == n) break;  // find path
    k++;
    int l = Diag(0, std::min(k, int(m)));
    int r = Diag(std::min(k, int(n)), 0);
    parlay::parallel_for(l, r + 1, [&](int id) {
      int t = -1;
      if (max_row[id] != -1) {
        int i = max_row[id];
        int j = i + m - id;
        if (i == n || j == m) {
          t = i;
        } else {
          int get_lcp =
              block_query_lcp(i + 1, j + 1, a, b, table_A, table_B, S_A, S_B);
          t = i + 1 + get_lcp;
        }
      }
      if (id > 0 && max_row[id - 1] != -1) {
        int i = max_row[id - 1];
        int j = i + m - id + 1;
        if (i == n) {
          t = n;
        } else {
          t = std::max(t, i + 1 +
                              block_query_lcp(i + 1, j, a, b, table_A, table_B,
                                              S_A, S_B));
        }
      }
      if (id < n + m && max_row[id + 1] != -1) {
        int i = max_row[id + 1];
        int j = i + m - id - 1;
        if (j == m) {
          t = std::max(t, i);
        } else {
          t = std::max(t, i + block_query_lcp(i, j + 1, a, b, table_A, table_B,
                                              S_A, S_B));
        }
      }
      // assert(t <= n);
      temp[id] = t;
    });
    parlay::parallel_for(l, r + 1,
                         [&](int id) { max_row[id] = std::min(temp[id], id); });
  }
  return k;
}
