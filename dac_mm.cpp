#include <parlay/delayed_sequence.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include <algorithm>

using namespace parlay;
using namespace std;

template <typename T, typename s_size_t = uint32_t>
class DAC_MM {
  static constexpr s_size_t MAX_VAL = std::numeric_limits<s_size_t>::max() / 2;
  struct Vector {
    sequence<s_size_t> &seq;
    const size_t c;
    Vector(sequence<s_size_t> &_seq, size_t _c) : seq(_seq), c(_c) {}
    s_size_t &operator[](size_t y) {
      if (y + c >= seq.size()) {
        printf("size: %zu, accessing %zu\n", seq.size(), y + c);
      }
      assert(y + c < seq.size());
      return seq[y + c];
    }
  };
  struct Matrix {
    sequence<sequence<s_size_t>> &seq;
    const size_t r;
    const size_t c;
    Matrix(sequence<sequence<s_size_t>> &_seq, size_t _r, size_t _c)
        : seq(_seq), r(_r), c(_c) {}
    Matrix(Matrix &_m, size_t _r, size_t _c)
        : seq(_m.seq), r(_m.r + _r), c(_m.c + _c) {}
    Vector operator[](size_t x) {
      if (x + r >= seq.size()) {
        printf("size: %zu, accessing %zu\n", seq.size(), x + r);
      }
      assert(x + r < seq.size());
      return Vector(seq[x + r], c);
    }

    Vector operator[](size_t x) const {
      assert(x + r < seq.size());
      return Vector(seq[x + r], c);
    }
  };

  const sequence<T> &A;
  const sequence<T> &B;

  DAC_MM(const sequence<T> &_A, const sequence<T> &_B) : A(_A), B(_B) {}

  size_t get_pow2(size_t x) {
    size_t ret = 1;
    while ((ret << 1) < x) {
      ret <<= 1;
    }
    return ret;
  }

  void merge_horizontal(Matrix left, Matrix right, Matrix ret, Matrix theta,
                        size_t n1, size_t n2, size_t k) {
    parallel_for(0, n1, [&](size_t i) {
      parallel_for(0, n1 - k, [&](size_t j) { ret[i][j] = left[i][j]; });
    });
    parallel_for(0, n2 - k, [&](size_t i) {
      parallel_for(0, n2, [&](size_t j) {
        ret[n1 - 1 + i][n1 - k - 1 + j] = right[k + i][j];
      });
    });
    parallel_for(0, n2 - k, [&](size_t i) {
      parallel_for(0, n1 - k, [&](size_t j) { ret[n1 - 1 + i][j] = MAX_VAL; });
    });

    auto compute = [&](size_t i, size_t j, size_t l, size_t r) {
      if (i + (n2 - k - 1) < j) {
        ret[i][n1 - k - 1 + j] = MAX_VAL;
        theta[i][j] = MAX_VAL;
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
      s_size_t mn_p = *min_element(perm, [&](s_size_t a, s_size_t b) {
        return get_val(a) < get_val(b);
      });
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
    theta[0][q] = theta[p][0] = MAX_VAL;
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
  }

  void merge_vertical(Matrix up, Matrix down, Matrix ret, Matrix theta,
                      size_t n1, size_t n2, size_t k) {
    parallel_for(0, n1, [&](size_t i) {
      parallel_for(0, n1 - k, [&](size_t j) {
        ret[n2 - k - 1 + i][n2 - 1 + j] = up[i][k + j];
      });
    });
    parallel_for(0, n2 - k, [&](size_t i) {
      parallel_for(0, n2, [&](size_t j) { ret[i][j] = down[i][j]; });
    });
    parallel_for(0, n2 - k, [&](size_t i) {
      parallel_for(0, n1 - k, [&](size_t j) { ret[i][n2 - 1 + j] = MAX_VAL; });
    });

    auto compute = [&](size_t i, size_t j, size_t l, size_t r) {
      if (i > j + (n1 - k - 1)) {
        ret[n2 - k - 1 + i][j] = MAX_VAL;
        theta[i][j] = MAX_VAL;
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
      s_size_t mn_p = *min_element(perm, [&](s_size_t a, s_size_t b) {
        return get_val(a) < get_val(b);
      });
      ret[n2 - k - 1 + i][j] = get_val(mn_p);
      theta[i][j] = mn_p;
    };
    auto compute_odd_even = [&](size_t p, size_t q) {
      for (size_t i = p; i < n1; i += p * 2) {
        for (size_t j = 0; j < n2; j += q * 2) {
          size_t s = theta[i - p][j];
          size_t t = i + p >= n1 ? k : theta[i + p][j];
          compute(i, j, s, t);
        }
      }
    };
    auto compute_even_odd = [&](size_t p, size_t q) {
      for (size_t i = 0; i < n1; i += p * 2) {
        for (size_t j = q; j < n2; j += q * 2) {
          size_t s = theta[i][j - q];
          size_t t = j + q >= n2 ? k : theta[i][j + q];
          compute(i, j, s, t);
        }
      }
    };
    auto compute_odd_odd = [&](size_t p, size_t q) {
      for (size_t i = p; i < n1; i += p * 2) {
        for (size_t j = q; j < n2; j += q * 2) {
          size_t s = theta[i][j - q];
          size_t t = j + q >= n2 ? k : theta[i][j + q];
          compute(i, j, s, t);
        }
      }
    };
    size_t p = get_pow2(n1), q = get_pow2(n2);
    compute(0, 0, 0, k);
    theta[0][q] = theta[p][0] = MAX_VAL;
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
  }

  void solve_r(size_t i, size_t n, size_t j, size_t m, Matrix dist,
               Matrix theta, Matrix tmp,
               const sequence<sequence<tuple<size_t, size_t, size_t, size_t,
                                             size_t, size_t>>> &dp) {
    size_t n1 = (n + 1) / 2, n2 = n - n1;
    size_t m1 = (m + 1) / 2, m2 = m - m1;
    if (n == 1 && m == 1) {
      dist[0][0] = dist[2][2] = 0;
      dist[0][1] = dist[1][0] = dist[1][2] = dist[2][1] = 1;
      dist[1][1] = (A[i] != B[j]);
    } else if (n == 1) {
      auto left = Matrix(dist, 0, 0);
      auto right = Matrix(dist, 0, get<1>(dp[n][m1]));
      par_do([&]() { solve_r(i, n, j, m1, left, theta, tmp, dp); },
             [&]() {
               solve_r(i, n, j + m1, m2, right,
                       Matrix(theta, 0, get<3>(dp[n][m1])),
                       Matrix(tmp, 0, get<5>(dp[n][m1])), dp);
             });
      merge_horizontal(left, right, Matrix(tmp, 0, 0), theta, n + m1 + 1,
                       n + m2 + 1, n);
      parallel_for(0, n + m + 1, [&](size_t x) {
        parallel_for(0, n + m + 1, [&](size_t y) { dist[x][y] = tmp[x][y]; });
      });
    } else if (m == 1) {
      auto up = Matrix(dist, 0, 0);
      auto down = Matrix(dist, get<0>(dp[n1][m]), 0);
      par_do([&]() { solve_r(i, n1, j, m, up, theta, tmp, dp); },
             [&]() {
               solve_r(i + n1, n2, j, m, down,
                       Matrix(theta, get<2>(dp[n1][m]), 0),
                       Matrix(tmp, get<4>(dp[n1][m]), 0), dp);
             });
      merge_vertical(up, down, Matrix(tmp, 0, 0), theta, n1 + m + 1, n2 + m + 1,
                     m);
      parallel_for(0, n + m + 1, [&](size_t x) {
        parallel_for(0, n + m + 1, [&](size_t y) { dist[x][y] = tmp[x][y]; });
      });
    } else {
      size_t len1, len2, len3, len4, len5, len6;
      std::tie(len1, len2, len3, len4, len5, len6) = dp[n1][m1];
      auto upper_left = Matrix(dist, 0, 0);
      auto bottom_left = Matrix(dist, len1, 0);
      auto upper_right = Matrix(dist, 0, len2);
      auto bottom_right = Matrix(dist, len1, len2);
      par_do(
          [&]() {
            par_do([&]() { solve_r(i, n1, j, m1, upper_left, theta, tmp, dp); },
                   [&]() {
                     solve_r(i + n1, n2, j, m1, bottom_left,
                             Matrix(theta, len3, 0), Matrix(tmp, len5, 0), dp);
                   });
          },
          [&]() {
            par_do(
                [&]() {
                  solve_r(i, n1, j + m1, m2, upper_right,
                          Matrix(theta, 0, len4), Matrix(tmp, 0, len6), dp);
                },
                [&]() {
                  solve_r(i + n1, n2, j + m1, m2, bottom_right,
                          Matrix(theta, len3, len4), Matrix(tmp, len5, len6),
                          dp);
                });
          });
      par_do(
          [&]() {
            merge_vertical(upper_left, bottom_left, Matrix(tmp, 0, 0), theta,
                           n1 + m1 + 1, n2 + m1 + 1, m1);
          },
          [&]() {
            merge_vertical(upper_right, bottom_right, Matrix(tmp, 0, len6),
                           Matrix(theta, 0, len4), n1 + m2 + 1, n2 + m2 + 1,
                           m2);
          });
      merge_horizontal(Matrix(tmp, 0, 0), Matrix(tmp, 0, len6), dist, theta,
                       n + m1 + 1, n + m2 + 1, n);
    }
  }

  size_t solve() {
    size_t n = A.size(), m = B.size();
    if (n == 0) {
      return m;
    } else if (m == 0) {
      return n;
    }

    // DP for matrix size
    auto dp = sequence<
        sequence<tuple<size_t, size_t, size_t, size_t, size_t, size_t>>>(
        n + 1,
        sequence<tuple<size_t, size_t, size_t, size_t, size_t, size_t>>(m + 1));
    for (size_t i = 1; i <= n; i++) {
      for (size_t j = 1; j <= m; j++) {
        if (i == 1 && j == 1) {
          dp[i][j] = make_tuple(3, 3, 3, 3, 4, 4);
        } else if (i == 1) {
          auto [n1, m1, u1, v1, p1, q1] = dp[i][(j + 1) / 2];
          auto [n2, m2, u2, v2, p2, q2] = dp[i][j / 2];
          dp[i][j] =
              make_tuple(max(n1, n2), m1 + m2, max(u1 + u2, (j + 1) / 2 + 2),
                         max(v1 + v2, j / 2 + 2), max(p1 + p2, i + j + 1),
                         max(q1 + q2, i + j + 1));
        } else if (j == 1) {
          auto [n1, m1, u1, v1, p1, q1] = dp[(i + 1) / 2][j];
          auto [n2, m2, u2, v2, p2, q2] = dp[i / 2][j];
          dp[i][j] =
              make_tuple(n1 + n2, max(m1, m2), max(u1 + u2, (i + 1) / 2 + 2),
                         max(v1 + v2, i / 2 + 2), max(p1 + p2, i + j + 1),
                         max(q1 + q2, i + j + 1));
        } else {
          auto [n1, m1, u1, v1, p1, q1] = dp[(i + 1) / 2][(j + 1) / 2];
          auto [n2, m2, u2, v2, p2, q2] = dp[i / 2][(j + 1) / 2];
          auto [n3, m3, u3, v3, p3, q3] = dp[(i + 1) / 2][j / 2];
          auto [n4, m4, u4, v4, p4, q4] = dp[i / 2][j / 2];
          auto p = max(n1, n3) + max(n2, n4);
          auto q = max(m1, m2) + max(m3, m4);
          auto w = max({max(u1, u3) + max(u2, u4),
                        (i + 1) / 2 + (j + 1) / 2 + 1, i + (j + 1) / 2 + 1});
          auto x = max(
              {max(v1, v2) + max(v3, v4), i / 2 * 2 + j + 2, i + j / 2 + 1});
          auto y = max(max(p1, p3) + max(p2, p4), i + (j + 1) / 2 + 1);
          auto z = max(max(q1, q2) + max(q3, q4), i + j / 2 + 1);
          dp[i][j] = make_tuple(p, q, w, x, y, z);
        }
        std::get<0>(dp[i][j]) = max(std::get<0>(dp[i][j]), i + j + 1);
        std::get<1>(dp[i][j]) = max(std::get<1>(dp[i][j]), i + j + 1);
      }
    }

    auto [N1, N2, N3, N4, N5, N6] = dp[n][m];
    auto dist =
        sequence<sequence<s_size_t>>(N1, sequence<s_size_t>(N2, MAX_VAL));
    auto theta =
        sequence<sequence<s_size_t>>(N3, sequence<s_size_t>::uninitialized(N4));
    auto tmp =
        sequence<sequence<s_size_t>>(N5, sequence<s_size_t>(N6, MAX_VAL));
    solve_r(0, n, 0, m, Matrix(dist, 0, 0), Matrix(theta, 0, 0),
            Matrix(tmp, 0, 0), dp);
    return dist[n][m];
  }
};

template class DAC_MM<uint32_t>;
