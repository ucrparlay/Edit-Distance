#ifndef hash_block_parallel_hpp
#define hash_block_parallel_hpp

#include "utils.h"
using namespace std;

// auxiliary function for power x^p
int mypower(int x, int p) {
  int res = 1;
  for (int i = 0; i < p; i++) {
    res *= x;
  }
  return res;
}

// build a table
template <typename T>
void build(const T seq, vector<vector<int>> &table_seq, size_t block_size) {
  size_t k = seq.size() / block_size;
  // pre-computed power table [p^(block_size), p^(2 * block_size), ... p^(logk *
  // block_size)]
  int LOG2_k = FASTLOG2(k);
  vector<int> block_power_table;
  for (int i = 0; i < LOG2_k; i++) {
    block_power_table.push_back(mypower(PRIME_BASE, int(block_size * (i + 1))));
  }
  table_seq.resize(LOG2_k + 1);
  // pre-computed log table [1, 2, 4, ..., logk]
  vector<int> block_log_table;
  for (int i = 0; i < LOG2_k + 1; i++) {
    block_log_table.push_back(1 << i);
  }

  // for the first dim of the table
  vector<int> aux_inside_block_table;
  for (int i = 0; i <= (int)(block_size); i++) {
    aux_inside_block_table.push_back(mypower(PRIME_BASE, i));
  }
  table_seq[0].resize(k);
  parlay::parallel_for(0, k, [&](int j) {
    // table_seq[0][j] = hash_value(seq, j * block_size, (j + 1) * block_size -
    // 1);
    table_seq[0][j] = 0;
    for (int b_j = j * block_size; b_j <= (int)((j + 1) * block_size - 1);
         b_j++) {
      table_seq[0][j] +=
          (aux_inside_block_table[(j + 1) * block_size - 1 - b_j] *
           (int)seq[b_j]);
    }
  });
  // for the remaining dims of the table
  for (int i = 1; i < LOG2_k + 1; i++) {
    table_seq[i].resize(k - block_log_table[i] + 1);
    parlay::parallel_for(0, k - block_log_table[i] + 1, [&](int j) {
      table_seq[i][j] = table_seq[i - 1][j] * block_power_table[i - 1] +
                        table_seq[i - 1][j + block_log_table[i - 1]];
    });
  }
}

// function to build the two tables. `n` is
// the length of the initial sequence size.
template <typename T>
size_t construct_table(T A, T B, vector<vector<int>> &table_A,
                       vector<vector<int>> &table_B, vector<int> &suffix_a,
                       vector<int> &suffix_b,
                       vector<int> &auxiliary_single_power_table, size_t n) {
  // logn
  int BLOCK_SIZE = FASTLOG2(n);
  // build the powertable
  if (BLOCK_SIZE == 0) {
    BLOCK_SIZE = 1;
  }
  build(A, table_A, BLOCK_SIZE);
  build(B, table_B, BLOCK_SIZE);

  // build the suffix
  // auxiliary power table
  vector<int> block_power_table;
  for (int i = 0; i < BLOCK_SIZE; i++) {
    block_power_table.push_back(mypower(PRIME_BASE, i));
  }

  size_t aux_size =
      std::max(table_A[0].size() * BLOCK_SIZE, table_B[0].size() * BLOCK_SIZE);
  auxiliary_single_power_table.resize(aux_size);
  auxiliary_single_power_table[0] = 1;
  for (int i = 1; i < aux_size; i++) {
    auxiliary_single_power_table[i] =
        PRIME_BASE * auxiliary_single_power_table[i - 1];
  };
  int a_actual_size = BLOCK_SIZE * int(A.size() / BLOCK_SIZE);
  int b_actual_size = BLOCK_SIZE * int(B.size() / BLOCK_SIZE);
  // prefix_a.resize(a_actual_size);
  // prefix_b.resize(b_actual_size);
  suffix_a.resize(a_actual_size);
  suffix_b.resize(b_actual_size);

  parlay::parallel_for(0, a_actual_size, [&](int i) {
    suffix_a[i] = 0;
    for (int k = i; k < (i / BLOCK_SIZE + 1) * int(BLOCK_SIZE); k++) {
      suffix_a[i] +=
          auxiliary_single_power_table[(i / BLOCK_SIZE + 1) * int(BLOCK_SIZE) -
                                       k - 1] *
          int(A[k]);
    }
  });

  parlay::parallel_for(0, b_actual_size, [&](int i) {
    suffix_b[i] = 0;
    for (int k = i; k < (i / BLOCK_SIZE + 1) * int(BLOCK_SIZE); k++) {
      suffix_b[i] +=
          auxiliary_single_power_table[(i / BLOCK_SIZE + 1) * int(BLOCK_SIZE) -
                                       k - 1] *
          int(B[k]);
    }
  });
  return BLOCK_SIZE;
}

bool compare_lcp(int p, int q, int z, vector<vector<int>> &table_A,
                 vector<vector<int>> &table_B, vector<int> &S_A,
                 vector<int> &S_B, vector<int> &aux_power_table, int t) {
  // size_t t = S_A.size() / table_A[0].size(); // block_size
  if (t == 0) {
    return false;
  }
  if ((p + (1 << z) * t) >= (int)(table_A[0].size() * t) ||
      (q + (1 << z) * t) >= (int)(table_B[0].size() * t)) {
    return false;
  }
  int next_block_A = (p / t) + 1;
  int next_block_B = (q / t) + 1;
  size_t rest_A_size = next_block_A * t - p;
  size_t rest_B_size = next_block_B * t - q;
  int hash_a_v;
  int hash_b_v;
  // A_hash_v = S_A[p] * p^{?} + Table_A[z][next_A] * p^{?} + P_A[nextA * t +
  // rest_B] need pre-computed prefix and suffix

  // if (p % t == 0) {
  //   hash_a_v = table_A[z][p / t];
  // } else if (z != 0) {
  //   hash_a_v = S_A[p] * aux_power_table[(1 << z) * t - rest_A_size] +
  //              (table_A[z][next_block_A]) -
  //              S_A[p + (1 << z) * t] / aux_power_table[rest_A_size];
  // } else {
  //   hash_a_v = S_A[p] * aux_power_table[t - rest_A_size] + P_A[p + t - 1];
  // }
  if (p % t == 0) {
    hash_a_v = table_A[z][p / t];
  } else {
    hash_a_v = S_A[p] * aux_power_table[(1 << z) * t - rest_A_size] +
               (table_A[z][next_block_A]) -
               S_A[p + (1 << z) * t] / aux_power_table[rest_A_size];
  }

  if (q % t == 0) {
    hash_b_v = table_B[z][q / t];
  } else {
    hash_b_v = S_B[q] * aux_power_table[(1 << z) * t - rest_B_size] +
               (table_B[z][next_block_B]) -
               S_B[q + (1 << z) * t] / aux_power_table[rest_B_size];
  }
  if (hash_a_v == hash_b_v) return true;
  return false;
}

// function for query the lcp from A[p] and B[q]
template <typename T>
int block_query_lcp(int p, int q, const T &A, const T &B,
                    vector<vector<int>> &table_A, vector<vector<int>> &table_B,
                    vector<int> &S_A, vector<int> &S_B,
                    vector<int> &aux_power_table, int t) {
  // find the possible block range point (omit the offset first)
  if (p >= (int)(A.size()) || q >= (int)(B.size())) {
    return 0;
  }

  if (A[p] != B[q]) {
    return 0;
  }
  int x = 0;
  // size_t t = S_A.size() / table_A[0].size(); // block_size
  while ((p + t * mypower(2, x) < (int)(table_A[0].size() * t)) &&
         (q + t * mypower(2, x) < (int)(table_B[0].size() * t))) {
    if (!compare_lcp(p, q, x, table_A, table_B, S_A, S_B, aux_power_table, t)) {
      break;
    }
    x++;
  }
  int pp = p;
  int qq = q;
  //    x -= 1; // now the x is the minimum LCP,
  // and the real lcp should locate in x and x + 1 block range
  if (x > 0) {
    pp = p + t * mypower(2, x - 1);
    qq = q + t * mypower(2, x - 1);
    int y = x - 1;
    while (y >= 0) {
      if (compare_lcp(pp, qq, y, table_A, table_B, S_A, S_B, aux_power_table,
                      t)) {
        pp += t * mypower(2, y);
        qq += t * mypower(2, y);
      }
      y -= 1;
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

  if (pp != p) {
    return (pp - p);
  } else {
    return 1;
  }
}

#endif