#ifndef hash_block_parallel_hpp
#define hash_block_parallel_hpp

#include "utils.h"
using namespace std;
constexpr int BLOCK_SIZE = 32;
constexpr int GRANULARITY = 512;

// auxiliary function for power x^p
// int mypower(int x, int p) {
// int res = 1;
// for (int i = 0; i < p; i++) {
// res *= x;
//}
// return res;
//}

int mypower(int a, int n) {
  int ret = 1;
  while (n) {
    if (n & 1) {
      ret = ret * a;
    }
    a = a * a;
    n >>= 1;
  }
  return ret;
}

// build a table
template <typename T>
void build(const T &seq, parlay::sequence<parlay::sequence<int>> &table_seq) {
  size_t k = seq.size() / BLOCK_SIZE;
  // pre-computed power table [p^(BLOCK_SIZE), p^(2 * BLOCK_SIZE), ... p^(logk *
  // BLOCK_SIZE)]
  int LOG2_k = FASTLOG2(k);
  table_seq = parlay::tabulate(k, [&](size_t i) {
    return parlay::sequence<int>::uninitialized(LOG2_k + 1);
  });
  // pre-computed log table [1, 2, 4, ..., logk]

  // for the first dim of the table
  parlay::parallel_for(0, k, [&](int j) {
    table_seq[j][0] = 0;
    for (int b_j = j * BLOCK_SIZE; b_j < (j + 1) * BLOCK_SIZE; b_j++) {
      table_seq[j][0] +=
          (mypower(PRIME_BASE, (j + 1) * BLOCK_SIZE - 1 - b_j) * seq[b_j]);
    }
  });
  // for the remaining dims of the table
  for (int i = 1; i < LOG2_k + 1; i++) {
    parlay::parallel_for(
        0, k - (1 << i) + 1,
        [&](int j) {
          table_seq[j][i] = table_seq[j][i - 1] * mypower(PRIME_BASE, i - 1) +
                            table_seq[j + (1 << (i - 1))][i - 1];
        },
        GRANULARITY);
  }
}

int qpow(int n) {
  int ret = 1;
  int a = PRIME_BASE;
  while (n) {
    if (n & 1) {
      ret = ret * a;
    }
    a = a * a;
    n >>= 1;
  }
  return ret;
}

// function to build the two tables. `n` is
// the length of the initial sequence size.
template <typename T>
void construct_table(T &A, T &B,
                     parlay::sequence<parlay::sequence<int>> &table_A,
                     parlay::sequence<parlay::sequence<int>> &table_B,
                     parlay::sequence<std::pair<int, int>> &pre_su_a,
                     parlay::sequence<std::pair<int, int>> &pre_su_b,
                     size_t n) {
  if (A.size() < BLOCK_SIZE || B.size() < BLOCK_SIZE) {
    return;
  }
  // logn

  /**
   * standard block size
   */
  parlay::internal::timer t;
  int BLOCK_SIZE_UPPER = FASTLOG2(n);

  /**
   * 32 / 64 ?
   */
  build(A, table_A);
  build(B, table_B);

  // build the suffix
  // auxiliary power table
  parlay::sequence<int> block_power_table;
  for (int i = 0; i < BLOCK_SIZE; i++) {
    block_power_table.push_back(mypower(PRIME_BASE, i));
  }

  int a_actual_size = BLOCK_SIZE * int(A.size() / BLOCK_SIZE);
  int b_actual_size = BLOCK_SIZE * int(B.size() / BLOCK_SIZE);
  // prefix_a.resize(a_actual_size);
  // prefix_b.resize(b_actual_size);

  pre_su_a =
      parlay::sequence<std::pair<int, int>>::uninitialized(a_actual_size);
  pre_su_b =
      parlay::sequence<std::pair<int, int>>::uninitialized(b_actual_size);
  // t.next("first");
  //  pre_su_a.resize(a_actual_size);
  //  pre_su_b.resize(b_actual_size);

  // for both prefix and suffix precomputing

  int num_blocks_a = a_actual_size / BLOCK_SIZE;
  int num_blocks_b = b_actual_size / BLOCK_SIZE;

  parlay::parallel_for(0, num_blocks_a, [&](size_t i) {
    int s = i * BLOCK_SIZE;
    int e = (i + 1) * BLOCK_SIZE;
    pre_su_a[s].first = A[s];
    for (int j = s + 1; j < e; j++) {
      pre_su_a[j].first = pre_su_a[j - 1].first * PRIME_BASE + A[j];
    }

    pre_su_a[e - 1].second = A[e - 1];
    for (int j = e - 2; j >= s; j--) {
      pre_su_a[j].second = pre_su_a[j + 1].second * PRIME_BASE + A[j];
    }
  });
  // t.next("second");

  parlay::parallel_for(0, num_blocks_b, [&](size_t i) {
    int s = i * BLOCK_SIZE;
    int e = (i + 1) * BLOCK_SIZE;
    pre_su_b[s].first = B[s];
    for (int j = s + 1; j < e; j++) {
      pre_su_b[j].first = pre_su_b[j - 1].first * PRIME_BASE + B[j];
    }

    pre_su_b[e - 1].second = B[e - 1];
    for (int j = e - 2; j >= s; j--) {
      pre_su_b[j].second = pre_su_b[j + 1].second * PRIME_BASE + B[j];
    }
  });
  // t.next("third");
}

template <typename T>
bool compare_lcp(int p, int q, int z,
                 const parlay::sequence<parlay::sequence<int>> &table_A,
                 const parlay::sequence<parlay::sequence<int>> &table_B,
                 const parlay::sequence<std::pair<int, int>> &S_A,
                 const parlay::sequence<std::pair<int, int>> &S_B, const T &A,
                 const T &B) {
  // size_t t = S_A.size() / table_A[0].size(); // BLOCK_SIZE
  constexpr int t = BLOCK_SIZE;
  if (t == 0) {
    return false;
  }
  int size_A = A.size() / t * t;
  int size_B = B.size() / t * t;
  if ((p + (1 << z) * t + t) >= size_A || (q + (1 << z) * t + t) >= size_B) {
    return false;
  }
  // if ((p + (1 << z) * t) >= (int)(A.size()) ||
  //     (q + (1 << z) * t) >= (int)(B.size())) {
  //   return false;
  // }
  int next_block_A = (p / t) + 1;
  int next_block_B = (q / t) + 1;
  size_t rest_A_size = next_block_A * t - p;
  size_t rest_B_size = next_block_B * t - q;
  int hash_a_v;
  int hash_b_v;

  if (p % t == 0) {
    hash_a_v = table_A[p / t][z] * qpow(t) + table_A[p / t + (1 << z)][0];
  } else {
    hash_a_v = S_A[p].second * qpow((1 << z) * t + t) +
               (table_A[next_block_A][z]) * qpow(t - rest_A_size) +
               S_A[p + (1 << z) * t + t].first;
  }

  if (q % t == 0) {
    hash_b_v = table_B[q / t][z] * qpow(t) + table_B[q / t + (1 << z)][0];
  } else {
    hash_b_v = S_B[q].second * qpow((1 << z) * t + t) +
               (table_B[next_block_B][z]) * qpow(t - rest_B_size) +
               S_B[q + (1 << z) * t + t].first;
  }
  if (hash_a_v == hash_b_v) return true;
  return false;
}

// function for query the lcp from A[p] and B[q]
template <typename T>
int block_query_lcp(int p, int q, const T &A, const T &B,
                    const parlay::sequence<parlay::sequence<int>> &table_A,
                    const parlay::sequence<parlay::sequence<int>> &table_B,
                    const parlay::sequence<std::pair<int, int>> &S_A,
                    const parlay::sequence<std::pair<int, int>> &S_B) {
  constexpr int t = BLOCK_SIZE;
  if (A.size() < t || B.size() < t) {
    int pp = p;
    int qq = q;
    while (pp < (int)A.size() && qq < (int)B.size() && A[pp] == B[qq]) {
      pp++;
      qq++;
    }
    return pp - p;
  }
  // find the possible block range point (omit the offset first)
  if (p >= (int)(A.size()) || q >= (int)(B.size())) {
    return 0;
  }

  if (A[p] != B[q]) {
    return 0;
  }
  int x = 0;
  // size_t t = S_A.size() / table_A[0].size(); // BLOCK_SIZE

  int size_A = A.size() / t * t;
  int size_B = B.size() / t * t;
  while ((p + (t << x) + t < size_A) && (q + (t << x) + t < size_B)) {
    if (!compare_lcp(p, q, x, table_A, table_B, S_A, S_B, A, B)) {
      break;
    }
    x++;
  }
  int pp = p;
  int qq = q;

  if (x > 0) {
    pp = p + (t << (x - 1)) + t;
    qq = q + (t << (x - 1)) + t;
    int y = x - 1;
    while (y >= 0) {
      if (compare_lcp(pp, qq, y, table_A, table_B, S_A, S_B, A, B)) {
        pp += (t << y) + t;
        qq += (t << y) + t;
      }
      y--;
    }
  } else {
    pp = p;
    qq = q;
  }

  while (pp < (int)(A.size()) && qq < (int)(B.size()) &&
         (int(A[pp]) == int(B[qq]))) {
    pp++;
    qq++;
  }

  return pp - p;
}

#endif
