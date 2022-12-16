#include <parlay/delayed_sequence.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include <algorithm>
#include <queue>

using namespace parlay;
using namespace std;

template <typename T, typename s_size_t = uint32_t>
class DAC_MM {
  static constexpr s_size_t MAX_VAL = std::numeric_limits<s_size_t>::max() / 2;
  static constexpr size_t BASE_CASE_SIZE = 64;
  template <typename It>
  struct Matrix {
    slice<It, It> seq;
    const size_t r;
    const size_t c;
    Matrix(slice<It, It> seq, size_t r, size_t c) : seq(seq), r(r), c(c) {}
    Matrix(Matrix &_m, size_t _r, size_t _c) : seq(_m.seq), r(_r), c(_c) {}
    slice<It, It> operator[](size_t x) { return seq.cut(x * c, x * (c + 1)); }
    slice<It, It> operator[](size_t x) const {
      return seq.cut(x * c, x * (c + 1));
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

  template <typename Iterator>
  void merge_horizontal(slice<Iterator, Iterator> _left,
                        slice<Iterator, Iterator> _right,
                        slice<Iterator, Iterator> _ret,
                        slice<Iterator, Iterator> _theta, size_t n1, size_t n2,
                        size_t k) {
    size_t n = n1 + n2 - k - 1;
    auto left = Matrix(_left, n1, n1);
    auto right = Matrix(_right, n2, n2);
    auto ret = Matrix(_ret, n, n);
    auto theta = Matrix(_theta, n1, n2);
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
      for (size_t i = 0; i < n2 - k; i++) {
        for (size_t j = 0; j < n1 - k; j++) {
          ret[n1 - 1 + i][j] = MAX_VAL;
        }
      }
      for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) {
          ret[i][n1 - k - 1 + j] = MAX_VAL;
          for (size_t o = 0; o < k + 1; o++) {
            if (left[i][n1 - k - 1 + o] != MAX_VAL && right[o][j] != MAX_VAL) {
              ret[i][n1 - k - 1 + j] =
                  min(ret[i][n1 - k - 1 + j],
                      left[i][n1 - k - 1 + o] + right[o][j]);
            }
          }
        }
      }
      return;
    }
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

  template <typename Iterator>
  void merge_vertical(slice<Iterator, Iterator> _up,
                      slice<Iterator, Iterator> _down,
                      slice<Iterator, Iterator> _ret,
                      slice<Iterator, Iterator> _theta, size_t n1, size_t n2,
                      size_t k) {
    size_t n = n1 + n2 - k - 1;
    auto up = Matrix(_up, n1, n1);
    auto down = Matrix(_down, n2, n2);
    auto ret = Matrix(_ret, n, n);
    auto theta = Matrix(_theta, n1, n2);
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
      for (size_t i = 0; i < n2 - k; i++) {
        for (size_t j = 0; j < n1 - k; j++) {
          ret[i][n2 - 1 + j] = MAX_VAL;
        }
      }
      for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) {
          ret[n2 - k - 1 + i][j] = MAX_VAL;
          for (size_t o = 0; o < k + 1; o++) {
            if (up[i][o] != MAX_VAL && down[n2 - k - 1 + o][j] != MAX_VAL) {
              ret[n2 - k - 1 + i][j] = min(ret[n2 - k - 1 + i][j],
                                           up[i][o] + down[n2 - k - 1 + o][j]);
            }
          }
        }
      }
      return;
    }
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

  size_t sqr(size_t x) { return x * x; }
  size_t get_size(size_t n, size_t m) { return sqr(2 * (n + m) - 1); }

  template <typename Iterator>
  s_size_t solve_r_serial(size_t i, size_t n, size_t j, size_t m,
                          slice<Iterator, Iterator> _dist,
                          slice<Iterator, Iterator> _theta,
                          slice<Iterator, Iterator> _tmp) {
    auto dist = Matrix<Iterator>(_dist, n + m + 1, n + m + 1);
    if (n + m + 1 < BASE_CASE_SIZE / 4) {
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
      return dist[n][m];
    }
    size_t n1 = (n + 1) / 2, n2 = n - n1;
    size_t m1 = (m + 1) / 2, m2 = m - m1;
    if (n == 1 && m == 1) {
      dist[0][0] = dist[2][2] = 0;
      dist[0][1] = dist[1][0] = dist[1][2] = dist[2][1] = 1;
      dist[1][1] = (A[i] != B[j]);
    } else if (n == 1) {
      size_t size1 = get_size(n, m1);
      size_t size2 = get_size(n, m2);
      size_t psize0 = 0;
      size_t psize1 = psize0 + size1;
      size_t psize2 = psize1 + size2;
      auto L = _dist.cut(psize0, psize1);
      auto R = _dist.cut(psize1, psize2);

      solve_r(i, n, j, m1, L, _theta.cut(0, size1), _tmp.cut(0, size1));
      solve_r(i, n, j + m1, m2, R, _theta.cut(size1, size1 + size2),
              _tmp.cut(size1, size1 + size2));

      merge_horizontal(L, R, _tmp, _theta, n + m1 + 1, n + m2 + 1, n);
      auto tmp = Matrix(_tmp, n + m + 1, n + m + 1);
      for (size_t x = 0; x < n + m + 1; x++) {
        for (size_t y = 0; y < n + m + 1; y++) {
          dist[x][y] = tmp[x][y];
        }
      }
    } else if (m == 1) {
      size_t size1 = get_size(n1, m);
      size_t size2 = get_size(n2, m);
      size_t psize0 = 0;
      size_t psize1 = psize0 + size1;
      size_t psize2 = psize1 + size2;
      auto U = _dist.cut(psize0, psize1);
      auto D = _dist.cut(psize1, psize2);

      solve_r(i, n1, j, m, U, _theta.cut(psize0, psize1),
              _tmp.cut(psize0, psize1));
      solve_r(i + n1, n2, j, m, D, _theta.cut(psize1, psize2),
              _tmp.cut(psize1, psize2));

      merge_vertical(U, D, _tmp, _theta, n1 + m + 1, n2 + m + 1, m);
      auto tmp = Matrix(_tmp, n + m + 1, n + m + 1);
      for (size_t x = 0; x < n + m + 1; x++) {
        for (size_t y = 0; y < n + m + 1; y++) {
          dist[x][y] = tmp[x][y];
        }
      }
    } else {
      size_t size1 = get_size(n1, m1);
      size_t size2 = get_size(n1, m2);
      size_t size3 = get_size(n2, m1);
      size_t size4 = get_size(n2, m2);
      size_t psize0 = 0;
      size_t psize1 = psize0 + size1;
      size_t psize2 = psize1 + size2;
      size_t psize3 = psize2 + size3;
      size_t psize4 = psize3 + size4;
      auto UL = _dist.cut(psize0, psize1);
      auto UR = _dist.cut(psize1, psize2);
      auto LL = _dist.cut(psize2, psize3);
      auto LR = _dist.cut(psize3, psize4);

      solve_r(i, n1, j, m1, UL, _theta.cut(psize0, psize1),
              _tmp.cut(psize0, psize1));
      solve_r(i, n1, j + m1, m2, UR, _theta.cut(psize1, psize2),
              _tmp.cut(psize1, psize2));
      solve_r(i + n1, n2, j, m1, LL, _theta.cut(psize2, psize3),
              _tmp.cut(psize2, size3));
      solve_r(i + n1, n2, j + m1, m2, LR, _theta.cut(psize3, psize4),
              _tmp.cut(psize3, psize4));

      auto theta1 = _theta.cut(psize0, psize2);
      auto theta2 = _theta.cut(psize2, psize4);
      auto theta3 = _theta.cut(psize0, psize4);
      auto tmp1 = _tmp.cut(psize0, psize2);
      auto tmp2 = _tmp.cut(psize2, psize4);

      merge_horizontal(UL, UR, tmp1, theta1, n1 + m1 + 1, n1 + m2 + 1, n1);
      merge_horizontal(LL, LR, tmp2, theta2, n2 + m1 + 1, n2 + m2 + 1, n2);
      merge_vertical(tmp1, tmp2, _dist, theta3, n1 + m + 1, n2 + m + 1, m);
    }
    return dist[n][m];
  }

  template <typename Iterator>
  s_size_t solve_r(size_t i, size_t n, size_t j, size_t m,
                   slice<Iterator, Iterator> _dist,
                   slice<Iterator, Iterator> _theta,
                   slice<Iterator, Iterator> _tmp) {
    if (n + m + 1 < BASE_CASE_SIZE) {
      return solve_r_serial(i, n, j, m, _dist, _theta, _tmp);
    }
    auto dist = Matrix<Iterator>(_dist, n + m + 1, n + m + 1);
    size_t n1 = (n + 1) / 2, n2 = n - n1;
    size_t m1 = (m + 1) / 2, m2 = m - m1;
    if (n == 1 && m == 1) {
      dist[0][0] = dist[2][2] = 0;
      dist[0][1] = dist[1][0] = dist[1][2] = dist[2][1] = 1;
      dist[1][1] = (A[i] != B[j]);
    } else if (n == 1) {
      size_t size1 = get_size(n, m1);
      size_t size2 = get_size(n, m2);
      size_t psize0 = 0;
      size_t psize1 = psize0 + size1;
      size_t psize2 = psize1 + size2;
      auto L = _dist.cut(psize0, psize1);
      auto R = _dist.cut(psize1, psize2);

      par_do(
          [&]() {
            solve_r(i, n, j, m1, L, _theta.cut(0, size1), _tmp.cut(0, size1));
          },
          [&]() {
            solve_r(i, n, j + m1, m2, R, _theta.cut(size1, size1 + size2),
                    _tmp.cut(size1, size1 + size2));
          });

      merge_horizontal(L, R, _tmp, _theta, n + m1 + 1, n + m2 + 1, n);
      auto tmp = Matrix(_tmp, n + m + 1, n + m + 1);
      parallel_for(0, n + m + 1, [&](size_t x) {
        parallel_for(0, n + m + 1, [&](size_t y) { dist[x][y] = tmp[x][y]; });
      });
    } else if (m == 1) {
      size_t size1 = get_size(n1, m);
      size_t size2 = get_size(n2, m);
      size_t psize0 = 0;
      size_t psize1 = psize0 + size1;
      size_t psize2 = psize1 + size2;
      auto U = _dist.cut(psize0, psize1);
      auto D = _dist.cut(psize1, psize2);

      par_do(
          [&]() {
            solve_r(i, n1, j, m, U, _theta.cut(psize0, psize1),
                    _tmp.cut(psize0, psize1));
          },
          [&]() {
            solve_r(i + n1, n2, j, m, D, _theta.cut(psize1, psize2),
                    _tmp.cut(psize1, psize2));
          });

      merge_vertical(U, D, _tmp, _theta, n1 + m + 1, n2 + m + 1, m);
      auto tmp = Matrix(_tmp, n + m + 1, n + m + 1);
      parallel_for(0, n + m + 1, [&](size_t x) {
        parallel_for(0, n + m + 1, [&](size_t y) { dist[x][y] = tmp[x][y]; });
      });
    } else {
      size_t size1 = get_size(n1, m1);
      size_t size2 = get_size(n1, m2);
      size_t size3 = get_size(n2, m1);
      size_t size4 = get_size(n2, m2);
      size_t psize0 = 0;
      size_t psize1 = psize0 + size1;
      size_t psize2 = psize1 + size2;
      size_t psize3 = psize2 + size3;
      size_t psize4 = psize3 + size4;
      auto UL = _dist.cut(psize0, psize1);
      auto UR = _dist.cut(psize1, psize2);
      auto LL = _dist.cut(psize2, psize3);
      auto LR = _dist.cut(psize3, psize4);

      par_do(
          [&]() {
            par_do(
                [&]() {
                  solve_r(i, n1, j, m1, UL, _theta.cut(psize0, psize1),
                          _tmp.cut(psize0, psize1));
                },
                [&]() {
                  solve_r(i, n1, j + m1, m2, UR, _theta.cut(psize1, psize2),
                          _tmp.cut(psize1, psize2));
                });
          },
          [&]() {
            par_do(
                [&]() {
                  solve_r(i + n1, n2, j, m1, LL, _theta.cut(psize2, psize3),
                          _tmp.cut(psize2, size3));
                },
                [&]() {
                  solve_r(i + n1, n2, j + m1, m2, LR,
                          _theta.cut(psize3, psize4), _tmp.cut(psize3, psize4));
                });
          });

      auto theta1 = _theta.cut(psize0, psize2);
      auto theta2 = _theta.cut(psize2, psize4);
      auto theta3 = _theta.cut(psize0, psize4);
      auto tmp1 = _tmp.cut(psize0, psize2);
      auto tmp2 = _tmp.cut(psize2, psize4);

      par_do(
          [&]() {
            merge_horizontal(UL, UR, tmp1, theta1, n1 + m1 + 1, n1 + m2 + 1,
                             n1);
          },
          [&]() {
            merge_horizontal(LL, LR, tmp2, theta2, n2 + m1 + 1, n2 + m2 + 1,
                             n2);
          });
      merge_vertical(tmp1, tmp2, _dist, theta3, n1 + m + 1, n2 + m + 1, m);
    }
    // printf("i: %zu, j: %zu, n: %zu, m: %zu\n", i, j, n, m);
    // for (size_t x = 0; x < n + m + 1; x++) {
    // for (size_t y = 0; y < n + m + 1; y++) {
    // printf("%10u ", dist[x][y]);
    //}
    // puts("");
    //}
    // puts("");
    return dist[n][m];
  }

  size_t solve() {
    size_t n = A.size(), m = B.size();
    if (n == 0) {
      return m;
    } else if (m == 0) {
      return n;
    }
    // printf("n: %zu, m: %zu, size: %zu\n", n, m, get_size(n, m));

    size_t size = get_size(n, m);
    auto dist = sequence<s_size_t>(size, MAX_VAL);
    auto theta = sequence<s_size_t>::uninitialized(size);
    auto tmp = sequence<s_size_t>(size, MAX_VAL);
    return solve_r(0, n, 0, m, make_slice(dist), make_slice(theta),
                   make_slice(tmp));
  }
};

template class DAC_MM<uint32_t>;
