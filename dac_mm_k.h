#ifndef DAC_MM_K_H
#define DAC_MM_K_H
#include <parlay/delayed_sequence.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include <algorithm>
#include <queue>

#include "dac_mm.h"

using namespace parlay;
using namespace std;

template <typename Seq, typename s_size_t = uint32_t>
class DAC_MM_K {
  static constexpr s_size_t MAX_VAL = std::numeric_limits<s_size_t>::max() / 2;
  const Seq &A;
  const Seq &B;
  const int max_recursion;

 public:
  DAC_MM_K(const Seq &_A, const Seq &_B)
      : A(_A), B(_B), max_recursion(log2_up(num_workers() * 10)) {}

  void print_matrix(const sequence<sequence<s_size_t>> &mat, string s = "") {
    printf("Matrix %s: ", s.c_str());
    if (mat.size() != 0) {
      printf("%zu by %zu\n", mat.size(), mat[0].size());
    } else {
      printf("0 by 0\n");
    }
    for (size_t i = 0; i < mat.size(); i++) {
      for (size_t j = 0; j < mat[0].size(); j++) {
        printf("%10u%c", mat[i][j], " \n"[j + 1 == mat[0].size()]);
      }
    }
  }

  sequence<sequence<s_size_t>> merge_horizontal(
      const sequence<sequence<s_size_t>> &left,
      const sequence<sequence<s_size_t>> &right, size_t k,
      bool do_parallel) {
    size_t n1 = left.size(), n2 = right.size();
    if (n1 == 0) {
      return right;
    }
    if (n2 == 0) {
      return left;
    }
    size_t n = n1 + n2 - k;
    size_t granularity = do_parallel ? 0 : std::numeric_limits<long>::max();
    sequence<sequence<s_size_t>> ret(n, sequence<s_size_t>(n, MAX_VAL));
    if (n < BASE_CASE_SIZE) {
      for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n1 - k; j++) {
          ret[i][j] = left[i][j];
        }
      }
      for (size_t i = 0; i < n2 - k; i++) {
        for (size_t j = 0; j < n2; j++) {
          ret[n1 + i][n1 - k + j] = right[k + i][j];
        }
      }
      for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) {
          for (size_t o = 0; o < k; o++) {
            if (left[i][n1 - k + o] != MAX_VAL && right[o][j] != MAX_VAL) {
              ret[i][n1 - k + j] =
                  min(ret[i][n1 - k + j], left[i][n1 - k + o] + right[o][j]);
            }
          }
        }
      }
      return ret;
    }
    parallel_for(
        0, n1,
        [&](size_t i) {
          parallel_for(
              0, n1 - k, [&](size_t j) { ret[i][j] = left[i][j]; },
              granularity);
        },
        granularity);
    parallel_for(
        0, n2 - k,
        [&](size_t i) {
          parallel_for(
              0, n2,
              [&](size_t j) { ret[n1 + i][n1 - k + j] = right[k + i][j]; },
              granularity);
        },
        granularity);

    sequence<sequence<s_size_t>> theta(n1, sequence<s_size_t>(n2, MAX_VAL));
    auto compute = [&](size_t i, size_t j, size_t l, size_t r) {
      if (l == MAX_VAL) {
        l = 0;
      }
      if (r == MAX_VAL) {
        r = k - 1;
      }
      auto get_val = [&](s_size_t o) {
        if (left[i][n1 - k + o] != MAX_VAL && right[o][j] != MAX_VAL) {
          return left[i][n1 - k + o] + right[o][j];
        }
        return MAX_VAL;
      };
      auto perm =
          delayed_seq<s_size_t>(r - l + 1, [&](size_t id) { return id + l; });
      s_size_t mn_p = *min_element(perm, [&](s_size_t a, s_size_t b) {
        return get_val(a) < get_val(b);
      });
      ret[i][n1 - k + j] = get_val(mn_p);
      theta[i][j] = mn_p;
    };
    auto compute_odd_even = [&](size_t p, size_t q) {
      parallel_for(
          0, (n1 - p - 1) / (p * 2) + 1,
          [&](size_t i) {
            parallel_for(
                0, (n2 - 1) / (q * 2) + 1,
                [&](size_t j) {
                  size_t x = p * (2 * i + 1);
                  size_t y = q * (2 * j);
                  size_t s = theta[x - p][y];
                  size_t t = x + p >= n1 ? k - 1 : theta[x + p][y];
                  compute(x, y, s, t);
                },
                granularity);
          },
          granularity);
    };
    auto compute_even_odd = [&](size_t p, size_t q) {
      parallel_for(
          0, (n1 - 1) / (p * 2) + 1,
          [&](size_t i) {
            parallel_for(
                0, (n2 - q - 1) / (q * 2) + 1,
                [&](size_t j) {
                  size_t x = p * (2 * i);
                  size_t y = q * (2 * j + 1);
                  size_t s = theta[x][y - q];
                  size_t t = y + q >= n2 ? k - 1 : theta[x][y + q];
                  compute(x, y, s, t);
                },
                granularity);
          },
          granularity);
    };
    auto compute_odd_odd = [&](size_t p, size_t q) {
      parallel_for(
          0, (n1 - p - 1) / (p * 2) + 1,
          [&](size_t i) {
            parallel_for(
                0, (n2 - q - 1) / (q * 2) + 1,
                [&](size_t j) {
                  size_t x = p * (2 * i + 1);
                  size_t y = q * (2 * j + 1);
                  size_t s = theta[x][y - q];
                  size_t t = y + q >= n2 ? k - 1 : theta[x][y + q];
                  compute(x, y, s, t);
                },
                granularity);
          },
          granularity);
    };
    size_t p = get_pow2(n1), q = get_pow2(n2);
    compute(0, 0, 0, k - 1);
    while (p && q) {
      compute_odd_even(p, q);
      compute_even_odd(p, q);
      compute_odd_odd(p, q);
      p >>= 1;
      q >>= 1;
    }
    while (p) {
      compute_odd_even(p, 1);
      compute_odd_odd(p, 1);
      p >>= 1;
    }
    while (q) {
      compute_even_odd(1, q);
      compute_odd_odd(1, q);
      q >>= 1;
    }
    return ret;
  }

  sequence<sequence<s_size_t>> merge_vertical(
      const sequence<sequence<s_size_t>> &up,
      const sequence<sequence<s_size_t>> &down, size_t k,
      bool do_parallel) {
    size_t n1 = up.size(), n2 = down.size();
    if (n1 == 0) {
      return down;
    }
    if (n2 == 0) {
      return up;
    }
    size_t n = n1 + n2 - k;
    size_t granularity = do_parallel ? 0 : std::numeric_limits<long>::max();
    sequence<sequence<s_size_t>> ret(n, sequence<s_size_t>(n, MAX_VAL));
    if (n < BASE_CASE_SIZE) {
      for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n1 - k; j++) {
          ret[n2 - k + i][n2 + j] = up[i][k + j];
        }
      }
      for (size_t i = 0; i < n2 - k; i++) {
        for (size_t j = 0; j < n2; j++) {
          ret[i][j] = down[i][j];
        }
      }
      for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) {
          for (size_t o = 0; o < k; o++) {
            if (up[i][o] != MAX_VAL && down[n2 - k + o][j] != MAX_VAL) {
              ret[n2 - k + i][j] =
                  min(ret[n2 - k + i][j], up[i][o] + down[n2 - k + o][j]);
            }
          }
        }
      }
      return ret;
    }
    parallel_for(
        0, n1,
        [&](size_t i) {
          parallel_for(
              0, n1 - k,
              [&](size_t j) { ret[n2 - k + i][n2 + j] = up[i][k + j]; },
              granularity);
        },
        granularity);
    parallel_for(
        0, n2 - k,
        [&](size_t i) {
          parallel_for(
              0, n2, [&](size_t j) { ret[i][j] = down[i][j]; }, granularity);
        },
        granularity);

    sequence<sequence<s_size_t>> theta(n1, sequence<s_size_t>(n2, MAX_VAL));
    auto compute = [&](size_t i, size_t j, size_t l, size_t r) {
      if (l == MAX_VAL) {
        l = 0;
      }
      if (r == MAX_VAL) {
        r = k - 1;
      }
      auto get_val = [&](s_size_t o) {
        if (up[i][o] != MAX_VAL && down[n2 - k + o][j] != MAX_VAL) {
          return up[i][o] + down[n2 - k + o][j];
        }
        return MAX_VAL;
      };
      auto perm =
          delayed_seq<s_size_t>(r - l + 1, [&](size_t id) { return id + l; });
      s_size_t mn_p = *min_element(perm, [&](s_size_t a, s_size_t b) {
        return get_val(a) < get_val(b);
      });
      ret[n2 - k + i][j] = get_val(mn_p);
      theta[i][j] = mn_p;
    };
    auto compute_odd_even = [&](size_t p, size_t q) {
      parallel_for(
          0, (n1 - p - 1) / (p * 2) + 1,
          [&](size_t i) {
            parallel_for(
                0, (n2 - 1) / (q * 2) + 1,
                [&](size_t j) {
                  size_t x = p * (2 * i + 1);
                  size_t y = q * (2 * j);
                  size_t s = theta[x - p][y];
                  size_t t = x + p >= n1 ? k - 1 : theta[x + p][y];
                  compute(x, y, s, t);
                },
                granularity);
          },
          granularity);
    };
    auto compute_even_odd = [&](size_t p, size_t q) {
      parallel_for(
          0, (n1 - 1) / (p * 2) + 1,
          [&](size_t i) {
            parallel_for(
                0, (n2 - q - 1) / (q * 2) + 1,
                [&](size_t j) {
                  size_t x = p * (2 * i);
                  size_t y = q * (2 * j + 1);
                  size_t s = theta[x][y - q];
                  size_t t = y + q >= n2 ? k - 1 : theta[x][y + q];
                  compute(x, y, s, t);
                },
                granularity);
          },
          granularity);
    };
    auto compute_odd_odd = [&](size_t p, size_t q) {
      parallel_for(
          0, (n1 - p - 1) / (p * 2) + 1,
          [&](size_t i) {
            parallel_for(
                0, (n2 - q - 1) / (q * 2) + 1,
                [&](size_t j) {
                  size_t x = p * (2 * i + 1);
                  size_t y = q * (2 * j + 1);
                  size_t s = theta[x][y - q];
                  size_t t = y + q >= n2 ? k - 1 : theta[x][y + q];
                  compute(x, y, s, t);
                },
                granularity);
          },
          granularity);
    };
    size_t p = get_pow2(n1), q = get_pow2(n2);
    compute(0, 0, 0, k - 1);
    while (p && q) {
      compute_odd_even(p, q);
      compute_even_odd(p, q);
      compute_odd_odd(p, q);
      p >>= 1;
      q >>= 1;
    }
    while (p) {
      compute_odd_even(p, 1);
      compute_odd_odd(p, 1);
      p >>= 1;
    }
    while (q) {
      compute_even_odd(1, q);
      compute_odd_odd(1, q);
      q >>= 1;
    }
    return ret;
  }

  sequence<sequence<s_size_t>> solve_t(size_t i, size_t n, size_t j, size_t m,
                                       bool upper_right, int recursion) {
    bool do_parallel =
        (n + m + 1 >= BASE_CASE_SIZE) && (recursion <= max_recursion);
    size_t granularity = do_parallel ? 0 : std::numeric_limits<long>::max();
    assert(n == m);
    if (n == 0) {
      return sequence<sequence<s_size_t>>();
    }
    auto tmp = DAC_MM<Seq, s_size_t>(A, B).solve_r(i, n, j, m, recursion);
    auto dist =
        sequence<sequence<s_size_t>>(n + 1, sequence<s_size_t>(n + 1, MAX_VAL));
    if (upper_right) {
      parallel_for(
          0, n + 1,
          [&](size_t i) {
            parallel_for(
                0, n + 1, [&](size_t j) { dist[i][j] = tmp[i + n][j + n]; },
                granularity);
          },
          granularity);
    } else {  // lower_left
      parallel_for(
          0, n + 1,
          [&](size_t i) {
            parallel_for(
                0, n + 1, [&](size_t j) { dist[i][j] = tmp[i][j]; },
                granularity);
          },
          granularity);
    }
    return dist;
  }

  sequence<sequence<s_size_t>> solve_r(size_t i, size_t n, size_t j, size_t m,
                                       size_t k1, size_t k2,
                                       int recursion = 0) {
    bool do_parallel =
        (n + m + 1 >= BASE_CASE_SIZE) && (recursion <= max_recursion);
    size_t granularity = do_parallel ? 0 : std::numeric_limits<long>::max();
    if (n / 2 <= k1 || m / 2 <= k2) {
      auto tmp = DAC_MM<Seq, s_size_t>(A, B).solve_r(i, n, j, m, recursion + 1);
      if (n >= k1 && m >= k2) {
        auto dist = sequence<sequence<s_size_t>>(
            k1 + k2 + 1, sequence<s_size_t>(k1 + k2 + 1, MAX_VAL));
        parallel_for(
            0, k1 + k2 + 1,
            [&](size_t i) {
              parallel_for(
                  0, k1 + k2 + 1,
                  [&](size_t j) {
                    dist[i][j] = tmp[i + (n - k1)][j + (n - k1)];
                  },
                  granularity);
            },
            granularity);
        return dist;
      } else if (n >= k1) {
        auto dist = sequence<sequence<s_size_t>>(
            k1 + m + 1, sequence<s_size_t>(k1 + m + 1, MAX_VAL));
        parallel_for(
            0, k1 + m + 1,
            [&](size_t i) {
              parallel_for(
                  0, k1 + m + 1,
                  [&](size_t j) {
                    dist[i][j] = tmp[i + (n - k1)][j + (n - k1)];
                  },
                  granularity);
            },
            granularity);
        return dist;
      } else if (m >= k2) {
        auto dist = sequence<sequence<s_size_t>>(
            n + k2 + 1, sequence<s_size_t>(n + k2 + 1, MAX_VAL));
        parallel_for(
            0, n + k2 + 1,
            [&](size_t i) {
              parallel_for(
                  0, n + k2 + 1, [&](size_t j) { dist[i][j] = tmp[i][j]; },
                  granularity);
            },
            granularity);
        return dist;
      } else {
        return tmp;
      }
    }
    if (n + m + 1 < BASE_CASE_SIZE / 4) {
      auto tmp = DAC_MM<Seq, s_size_t>(A, B).solve_r_serial(i, n, j, m);
      auto dist = sequence<sequence<s_size_t>>(
          k1 + k2 + 1, sequence<s_size_t>(k1 + k2 + 1, MAX_VAL));
      for (size_t i = 0; i < k1 + k2 + 1; i++) {
        for (size_t j = 0; j < k1 + k2 + 1; j++) {
          dist[i][j] = tmp[i + (n - k1)][j + (n - k1)];
        }
      }
      return dist;
    }
    size_t n1 = n / 2, n2 = n - n1;
    size_t m1 = m / 2, m2 = m - m1;

    size_t s1 = n1 + k2 - m1;
    size_t s2 = m1 + k1 - n1;

    sequence<sequence<s_size_t>> UL, UR, LL, LR, U, D;
    par_do_if(
        do_parallel,
        [&]() {
          par_do_if(
              do_parallel,
              [&]() { UL = solve_r(i, n1, j, m1, k1, k2, recursion + 1); },
              [&]() {
                UR = solve_t(i + n1 - s1, s1, j + m1, s1, false, recursion + 1);
              });
        },
        [&]() {
          par_do_if(
              do_parallel,
              [&]() {
                LL = solve_t(i + n1, s2, j + m1 - s2, s2, true, recursion + 1);
              },
              [&]() {
                LR = solve_r(i + n1, n2, j + m1, m2, s2, s1, recursion + 1);
              });
        });
    par_do_if(
        do_parallel,
        [&]() { U = merge_horizontal(UL, UR, min(s1, n1) + 1, do_parallel); },
        [&]() { D = merge_horizontal(LL, LR, min(s2, n2) + 1, do_parallel); });
    auto dist =
        merge_vertical(U, D, min(s1, m2) + min(s2, m1) + 1, do_parallel);
    return dist;
  }

  size_t solve_k(size_t k) {
    size_t n = A.size(), m = B.size();
    if (k >= min(n, m)) {
      return DAC_MM<Seq, s_size_t>(A, B).solve();
    } else {
      auto dist = solve_r(0, n, 0, m, k, k);
      return dist[k][m + k - n];
    }
  }

  size_t solve() {
    size_t n = A.size(), m = B.size();
    if (n == 0) {
      return m;
    } else if (m == 0) {
      return n;
    }

    size_t l = max((size_t)1, (n > m) ? (n - m) : (m - n));
    size_t r = max(n, m);
    size_t ret = 0;
    while (l <= r) {
      ret = solve_k(l);
      if (ret <= l) {
        break;
      }
      l = min(l * 2, r);
    }
    return ret;
  }
};

#endif  // DAC_MM_K_H
