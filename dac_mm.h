#ifndef DAC_MM_H
#define DAC_MM_H
#include <parlay/delayed_sequence.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include <algorithm>
#include <queue>

using namespace parlay;
using namespace std;
static constexpr size_t BASE_CASE_SIZE = 0;

size_t get_pow2(size_t x) {
  size_t ret = 1;
  while ((ret << 1) < x) {
    ret <<= 1;
  }
  return ret;
}

template <typename s_size_t>
sequence<sequence<s_size_t>> merge_horizontal(
    const sequence<sequence<s_size_t>> &left,
    const sequence<sequence<s_size_t>> &right, size_t k) {
  size_t n1 = left.size(), n2 = right.size();
  if (n1 == 0) {
    return right;
  }
  if (n2 == 0) {
    return left;
  }
  size_t n = n1 + n2 - k - 1;
  static constexpr s_size_t MAX_VAL = std::numeric_limits<s_size_t>::max() / 2;
  sequence<sequence<s_size_t>> ret(n, sequence<s_size_t>(n, MAX_VAL));
  if (n < BASE_CASE_SIZE) {
    for (size_t i = 0; i < n1; i++) {
      for (size_t j = 0; j < n1 - k; j++) {
        ret[i][j] = left[i][j];
      }
    }
    for (size_t i = 0; i < n2 - k; i++) {
      for (size_t j = 0; j < n2; j++) {
        ret[n1 - 1 + i][n1 - k - 1 + j] = right[k + i][j];
      }
    }
    for (size_t i = 0; i < n1; i++) {
      for (size_t j = 0; j < n2; j++) {
        for (size_t o = 0; o < k + 1; o++) {
          if (left[i][n1 - k - 1 + o] != MAX_VAL && right[o][j] != MAX_VAL) {
            ret[i][n1 - k - 1 + j] = min(ret[i][n1 - k - 1 + j],
                                         left[i][n1 - k - 1 + o] + right[o][j]);
          }
        }
      }
    }
    return ret;
  }
  parallel_for(0, n1, [&](size_t i) {
    parallel_for(0, n1 - k, [&](size_t j) { ret[i][j] = left[i][j]; });
  });
  parallel_for(0, n2 - k, [&](size_t i) {
    parallel_for(0, n2, [&](size_t j) {
      ret[n1 - 1 + i][n1 - k - 1 + j] = right[k + i][j];
    });
  });

  sequence<sequence<s_size_t>> theta(n1, sequence<s_size_t>(n2, MAX_VAL));
  auto compute = [&](size_t i, size_t j, size_t l, size_t r) {
    if (i + (n2 - k - 1) < j) {
      return;
    }
    if (l == MAX_VAL) {
      l = 0;
    }
    if (r == MAX_VAL) {
      r = k;
    }
    auto get_val = [&](s_size_t o) {
      if (left[i][n1 - k - 1 + o] != MAX_VAL && right[o][j] != MAX_VAL) {
        return left[i][n1 - k - 1 + o] + right[o][j];
      }
      return MAX_VAL;
    };
    auto perm =
        delayed_seq<s_size_t>(r - l + 1, [&](size_t id) { return id + l; });
    s_size_t mn_p = *min_element(
        perm, [&](s_size_t a, s_size_t b) { return get_val(a) < get_val(b); });
    ret[i][n1 - k - 1 + j] = get_val(mn_p);
    theta[i][j] = mn_p;
  };
  auto compute_odd_even = [&](size_t p, size_t q) {
    parallel_for(0, (n1 - p - 1) / (p * 2) + 1, [&](size_t i) {
      parallel_for(0, (n2 - 1) / (q * 2) + 1, [&](size_t j) {
        size_t x = p * (2 * i + 1);
        size_t y = q * (2 * j);
        size_t s = theta[x - p][y];
        size_t t = x + p >= n1 ? k : theta[x + p][y];
        compute(x, y, s, t);
      });
    });
  };
  auto compute_even_odd = [&](size_t p, size_t q) {
    parallel_for(0, (n1 - 1) / (p * 2) + 1, [&](size_t i) {
      parallel_for(0, (n2 - q - 1) / (q * 2) + 1, [&](size_t j) {
        size_t x = p * (2 * i);
        size_t y = q * (2 * j + 1);
        size_t s = theta[x][y - q];
        size_t t = y + q >= n2 ? k : theta[x][y + q];
        compute(x, y, s, t);
      });
    });
  };
  auto compute_odd_odd = [&](size_t p, size_t q) {
    parallel_for(0, (n1 - p - 1) / (p * 2) + 1, [&](size_t i) {
      parallel_for(0, (n2 - q - 1) / (q * 2) + 1, [&](size_t j) {
        size_t x = p * (2 * i + 1);
        size_t y = q * (2 * j + 1);
        size_t s = theta[x][y - q];
        size_t t = y + q >= n2 ? k : theta[x][y + q];
        compute(x, y, s, t);
      });
    });
  };
  size_t p = get_pow2(n1), q = get_pow2(n2);
  compute(0, 0, 0, k);
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

template <typename s_size_t>
sequence<sequence<s_size_t>> merge_vertical(
    const sequence<sequence<s_size_t>> &up,
    const sequence<sequence<s_size_t>> &down, size_t k) {
  size_t n1 = up.size(), n2 = down.size();
  if (n1 == 0) {
    return down;
  }
  if (n2 == 0) {
    return up;
  }
  size_t n = n1 + n2 - k - 1;
  static constexpr s_size_t MAX_VAL = std::numeric_limits<s_size_t>::max() / 2;
  sequence<sequence<s_size_t>> ret(n, sequence<s_size_t>(n, MAX_VAL));
  if (n < BASE_CASE_SIZE) {
    for (size_t i = 0; i < n1; i++) {
      for (size_t j = 0; j < n1 - k; j++) {
        ret[n2 - k - 1 + i][n2 - 1 + j] = up[i][k + j];
      }
    }
    for (size_t i = 0; i < n2 - k; i++) {
      for (size_t j = 0; j < n2; j++) {
        ret[i][j] = down[i][j];
      }
    }
    for (size_t i = 0; i < n1; i++) {
      for (size_t j = 0; j < n2; j++) {
        for (size_t o = 0; o < k + 1; o++) {
          if (up[i][o] != MAX_VAL && down[n2 - k - 1 + o][j] != MAX_VAL) {
            ret[n2 - k - 1 + i][j] =
                min(ret[n2 - k - 1 + i][j], up[i][o] + down[n2 - k - 1 + o][j]);
          }
        }
      }
    }
    return ret;
  }
  parallel_for(0, n1, [&](size_t i) {
    parallel_for(0, n1 - k, [&](size_t j) {
      ret[n2 - k - 1 + i][n2 - 1 + j] = up[i][k + j];
    });
  });
  parallel_for(0, n2 - k, [&](size_t i) {
    parallel_for(0, n2, [&](size_t j) { ret[i][j] = down[i][j]; });
  });

  sequence<sequence<s_size_t>> theta(n1, sequence<s_size_t>(n2, MAX_VAL));
  auto compute = [&](size_t i, size_t j, size_t l, size_t r) {
    if (i > j + (n1 - k - 1)) {
      return;
    }
    if (l == MAX_VAL) {
      l = 0;
    }
    if (r == MAX_VAL) {
      r = k;
    }
    auto get_val = [&](s_size_t o) {
      if (up[i][o] != MAX_VAL && down[n2 - k - 1 + o][j] != MAX_VAL) {
        return up[i][o] + down[n2 - k - 1 + o][j];
      }
      return MAX_VAL;
    };
    auto perm =
        delayed_seq<s_size_t>(r - l + 1, [&](size_t id) { return id + l; });
    s_size_t mn_p = *min_element(
        perm, [&](s_size_t a, s_size_t b) { return get_val(a) < get_val(b); });
    ret[n2 - k - 1 + i][j] = get_val(mn_p);
    theta[i][j] = mn_p;
  };
  auto compute_odd_even = [&](size_t p, size_t q) {
    parallel_for(0, (n1 - p - 1) / (p * 2) + 1, [&](size_t i) {
      parallel_for(0, (n2 - 1) / (q * 2) + 1, [&](size_t j) {
        size_t x = p * (2 * i + 1);
        size_t y = q * (2 * j);
        size_t s = theta[x - p][y];
        size_t t = x + p >= n1 ? k : theta[x + p][y];
        compute(x, y, s, t);
      });
    });
  };
  auto compute_even_odd = [&](size_t p, size_t q) {
    parallel_for(0, (n1 - 1) / (p * 2) + 1, [&](size_t i) {
      parallel_for(0, (n2 - q - 1) / (q * 2) + 1, [&](size_t j) {
        size_t x = p * (2 * i);
        size_t y = q * (2 * j + 1);
        size_t s = theta[x][y - q];
        size_t t = y + q >= n2 ? k : theta[x][y + q];
        compute(x, y, s, t);
      });
    });
  };
  auto compute_odd_odd = [&](size_t p, size_t q) {
    parallel_for(0, (n1 - p - 1) / (p * 2) + 1, [&](size_t i) {
      parallel_for(0, (n2 - q - 1) / (q * 2) + 1, [&](size_t j) {
        size_t x = p * (2 * i + 1);
        size_t y = q * (2 * j + 1);
        size_t s = theta[x][y - q];
        size_t t = y + q >= n2 ? k : theta[x][y + q];
        compute(x, y, s, t);
      });
    });
  };
  size_t p = get_pow2(n1), q = get_pow2(n2);
  compute(0, 0, 0, k);
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

template <typename Seq, typename s_size_t = uint32_t>
class DAC_MM {
  const Seq &A;
  const Seq &B;

 public:
  DAC_MM(const Seq &_A, const Seq &_B) : A(_A), B(_B) {}

  sequence<sequence<s_size_t>> solve_r_serial(size_t i, size_t n, size_t j,
                                              size_t m) {
    auto dist =
        sequence<sequence<s_size_t>>(n + m + 1, sequence<s_size_t>(n + m + 1));
    static constexpr s_size_t MAX_VAL =
        std::numeric_limits<s_size_t>::max() / 2;
    for (size_t x = 0; x < n + m + 1; x++) {
      s_size_t d[n + 1][m + 1];
      for (size_t u = 0; u < n + 1; u++) {
        for (size_t v = 0; v < m + 1; v++) {
          d[u][v] = MAX_VAL;
        }
      }
      s_size_t s = (x <= n ? n - x : 0);
      s_size_t t = (x <= n ? 0 : x - n);
      queue<pair<s_size_t, s_size_t>> q;
      q.push(make_pair(s, t));
      d[s][t] = 0;
      while (!q.empty()) {
        auto [u, v] = q.front();
        q.pop();
        if (u + 1 <= n && v + 1 <= m) {
          s_size_t w = A[i + u] == B[j + v] ? 0 : 1;
          if (d[u + 1][v + 1] > d[u][v] + w) {
            d[u + 1][v + 1] = d[u][v] + w;
            q.push(make_pair(u + 1, v + 1));
          }
        }
        if (u + 1 <= n && d[u + 1][v] > d[u][v] + 1) {
          d[u + 1][v] = d[u][v] + 1;
          q.push(make_pair(u + 1, v));
        }
        if (v + 1 <= m && d[u][v + 1] > d[u][v] + 1) {
          d[u][v + 1] = d[u][v] + 1;
          q.push(make_pair(u, v + 1));
        }
      }
      for (size_t y = 0; y < n + m + 1; y++) {
        s_size_t s = (y <= m ? n : n - (y - m));
        s_size_t t = (y <= m ? y : m);
        dist[x][y] = d[s][t];
      }
    }
    return dist;
  }

  sequence<sequence<s_size_t>> solve_r(size_t i, size_t n, size_t j, size_t m) {
    if (n == 0 || m == 0) {
      return sequence<sequence<s_size_t>>();
    }
    if (n + m + 1 < BASE_CASE_SIZE / 4) {
      return solve_r_serial(i, n, j, m);
    }
    bool do_parallel = (n + m + 1 >= BASE_CASE_SIZE);
    size_t n1 = (n + 1) / 2, n2 = n - n1;
    size_t m1 = (m + 1) / 2, m2 = m - m1;
    static constexpr s_size_t MAX_VAL =
        std::numeric_limits<s_size_t>::max() / 2;
    if (n == 1 && m == 1) {
      auto dist = sequence<sequence<s_size_t>>(
          n + m + 1, sequence<s_size_t>(n + m + 1, MAX_VAL));
      dist[0][0] = dist[2][2] = 0;
      dist[0][1] = dist[1][0] = dist[1][2] = dist[2][1] = 1;
      dist[1][1] = (A[i] != B[j]);
      return dist;
    } else if (m == 1) {
      sequence<sequence<s_size_t>> U, D;
      par_do_if(
          do_parallel, [&]() { U = solve_r(i, n1, j, m); },
          [&]() { D = solve_r(i + n1, n2, j, m); });

      auto dist = merge_vertical(U, D, m);
      return dist;
    } else {
      sequence<sequence<s_size_t>> UL, UR, LL, LR, U, D;
      par_do_if(
          do_parallel,
          [&]() {
            par_do_if(
                do_parallel, [&]() { UL = solve_r(i, n1, j, m1); },
                [&]() { UR = solve_r(i, n1, j + m1, m2); });
          },
          [&]() {
            par_do_if(
                do_parallel, [&]() { LL = solve_r(i + n1, n2, j, m1); },
                [&]() { LR = solve_r(i + n1, n2, j + m1, m2); });
          });

      par_do_if(
          do_parallel, [&]() { U = merge_horizontal(UL, UR, n1); },
          [&]() { D = merge_horizontal(LL, LR, n2); });
      auto dist = merge_vertical(U, D, m);
      return dist;
    }
  }

  size_t solve() {
    size_t n = A.size(), m = B.size();
    if (n == 0) {
      return m;
    } else if (m == 0) {
      return n;
    }
    auto dist = solve_r(0, n, 0, m);
    return dist[n][m];
  }
};

#endif  // DAC_MM_H
