#include <parlay/delayed_sequence.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

using namespace parlay;

template <typename T, typename s_size_t = uint32_t>
class DAC_MM {
  static constexpr s_size_t MAX_VAL = std::numeric_limits<s_size_t>::max() / 2;
  const sequence<T> &A;
  const sequence<T> &B;

  DAC_MM(const sequence<T> &_A, const sequence<T> &_B) : A(_A), B(_B) {}

  sequence<sequence<s_size_t>> merge_horizontal(
      sequence<sequence<s_size_t>> &left, sequence<sequence<s_size_t>> &right,
      size_t k) {
    size_t n1 = left.size(), n2 = right.size();
    size_t n = n1 + n2 - k - 1;
    sequence<sequence<s_size_t>> ret(n, sequence<s_size_t>(n, MAX_VAL));
    parallel_for(0, n1, [&](size_t i) {
      parallel_for(0, n1 - k, [&](size_t j) { ret[i][j] = left[i][j]; });
    });
    parallel_for(0, n2 - k, [&](size_t i) {
      parallel_for(0, n2, [&](size_t j) {
        ret[n1 - 1 + i][n1 - k - 1 + j] = right[k + i][j];
      });
    });

    parallel_for(0, n1, [&](size_t i) {
      parallel_for(0, n2, [&](size_t j) {
        ret[i][n1 - k - 1 + j] = reduce(
            delayed_seq<s_size_t>(k + 1,
                                  [&](size_t o) {
                                    return std::min(
                                        MAX_VAL,
                                        left[i][n1 - k - 1 + o] + right[o][j]);
                                  }),
            minm<s_size_t>());
      });
    });
    return ret;
  }

  sequence<sequence<s_size_t>> merge_vertical(
      sequence<sequence<s_size_t>> &up, sequence<sequence<s_size_t>> &down,
      size_t k) {
    size_t n1 = up.size(), n2 = down.size();
    size_t n = n1 + n2 - k - 1;
    sequence<sequence<s_size_t>> ret(n, sequence<s_size_t>(n, MAX_VAL));
    parallel_for(0, n1, [&](size_t i) {
      parallel_for(0, n1 - k, [&](size_t j) {
        ret[n2 - k - 1 + i][n2 - 1 + j] = up[i][k + j];
      });
    });
    parallel_for(0, n2 - k, [&](size_t i) {
      parallel_for(0, n2, [&](size_t j) { ret[i][j] = down[i][j]; });
    });
    parallel_for(0, n1, [&](size_t i) {
      parallel_for(0, n2, [&](size_t j) {
        ret[n2 - k - 1 + i][j] = reduce(
            delayed_seq<s_size_t>(
                k + 1,
                [&](size_t o) {
                  return std::min(MAX_VAL, up[i][o] + down[n2 - k - 1 + o][j]);
                }),
            minm<s_size_t>());
      });
    });
    return ret;
  }

  void print_matrix(size_t i, size_t n, size_t j, size_t m,
                    sequence<sequence<s_size_t>> &ret) {
    printf("i: %zu, n: %zu, j: %zu, m: %zu\n", i, n, j, m);
    for (size_t o = 0; o < ret.size(); o++) {
      for (size_t oo = 0; oo < ret[o].size(); oo++) {
        printf("%10u%c", ret[o][oo], " \n"[oo + 1 == ret[o].size()]);
      }
    }
    puts("");
  }

  sequence<sequence<s_size_t>> solve_r(size_t i, size_t n, size_t j, size_t m) {
    // printf("n: %zu, m: %zu\n", n, m);
    size_t n1 = n / 2, n2 = n - n1;
    size_t m1 = m / 2, m2 = m - m1;
    if (n == 1 && m == 1) {
      auto ret =
          sequence<sequence<s_size_t>>(3, sequence<s_size_t>(3, MAX_VAL));
      ret[0][0] = ret[2][2] = 0;
      ret[0][1] = ret[1][0] = ret[1][2] = ret[2][1] = 1;
      ret[1][1] = (A[i] != B[j]);
      // print_matrix(i, n, j, m, ret);
      return ret;
    } else if (n == 1) {
      sequence<sequence<s_size_t>> v1, v2;
      par_do([&]() { v1 = solve_r(i, n, j, m1); },
             [&]() { v2 = solve_r(i, n, j + m1, m2); });
      auto ret = merge_horizontal(v1, v2, n);
      // print_matrix(i, n, j, m, ret);
      return ret;
    } else if (m == 1) {
      sequence<sequence<s_size_t>> v1, v2;
      par_do([&]() { v1 = solve_r(i, n1, j, m); },
             [&]() { v2 = solve_r(i + n1, n2, j, m); });
      auto ret = merge_vertical(v1, v2, m);
      // print_matrix(i, n, j, m, ret);
      return ret;
    } else {
      sequence<sequence<s_size_t>> t1, t2, t3, t4;
      par_do(
          [&]() {
            par_do([&]() { t1 = solve_r(i, n1, j, m1); },
                   [&]() { t2 = solve_r(i + n1, n2, j, m1); });
          },
          [&]() {
            par_do([&]() { t3 = solve_r(i, n1, j + m1, m2); },
                   [&]() { t4 = solve_r(i + n1, n2, j + m1, m2); });
          });

      sequence<sequence<s_size_t>> v1, v2;
      par_do([&]() { v1 = merge_vertical(t1, t2, m1); },
             [&]() { v2 = merge_vertical(t3, t4, m2); });
      // print_matrix(i, n, j, m1, v1);
      // print_matrix(i, n, j + m1, m2, v2);
      auto ret = merge_horizontal(v1, v2, n);
      // print_matrix(i, n, j, m, ret);
      return ret;
    }
  }

  size_t solve() {
    size_t n = A.size(), m = B.size();
    auto distance = solve_r(0, n, 0, m);
    return distance[n][m];
  }
};

template class DAC_MM<uint32_t>;
