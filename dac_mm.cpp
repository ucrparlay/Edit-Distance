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

  sequence<sequence<s_size_t>> dist;
  sequence<sequence<s_size_t>> g_theta;

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
    // for (size_t i = 0; i < n1 + n2 - k - 1; i++) {
    // for (size_t j = 0; j < n1 + n2 - k - 1; j++) {
    // ret[i][j] = MAX_VAL;
    //}
    //}
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

    auto compute = [&](size_t i, size_t j, size_t l, size_t r) {
      if (i + k < j) {
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
      s_size_t mn_v = MAX_VAL;
      s_size_t mn_p = MAX_VAL;
      for (size_t o = l; o <= r; o++) {
        s_size_t sum = left[i][n1 - k - 1 + o] + right[o][j];
        if (sum < mn_v) {
          mn_v = sum;
          mn_p = o;
        }
      }
      ret[i][n1 - k - 1 + j] = mn_v;
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

  void merge_vertical(Matrix up, Matrix down, Matrix ret, Matrix theta,
                      size_t n1, size_t n2, size_t k) {
    // for (size_t i = 0; i < n1 + n2 - k - 1; i++) {
    // for (size_t j = 0; j < n1 + n2 - k - 1; j++) {
    // ret[i][j] = MAX_VAL;
    //}
    //}
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
    for (size_t i = 0; i < n2 - k - 1; i++) {
      for (size_t j = 0; j < n1 - k; j++) {
        ret[i][n2 - 1 + j] = MAX_VAL;
      }
    }

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
      s_size_t mn_v = MAX_VAL;
      s_size_t mn_p = MAX_VAL;
      for (size_t o = l; o <= r; o++) {
        s_size_t sum = up[i][o] + down[n2 - k - 1 + o][j];
        if (sum < mn_v) {
          mn_v = sum;
          mn_p = o;
        }
      }
      ret[n2 - k - 1 + i][j] = mn_v;
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

  void solve_r(size_t i, size_t n, size_t j, size_t m, size_t r, size_t c) {
    size_t n1 = (n + 1) / 2, n2 = n - n1;
    size_t m1 = (m + 1) / 2, m2 = m - m1;
    Matrix ret(dist, r, c);
    Matrix theta(g_theta, r, c);
    if (n == 1 && m == 1) {
      ret[0][0] = ret[2][2] = 0;
      ret[0][1] = ret[1][0] = ret[1][2] = ret[2][1] = 1;
      ret[1][1] = (A[i] != B[j]);
    } else if (n == 1) {
      auto left = Matrix(ret, 0, 0);
      auto right = Matrix(ret, 0, m1 * 3);
      solve_r(i, n, j, m1, r, c);
      solve_r(i, n, j + m1, m2, r, c + m1 * 3);
      sequence<sequence<s_size_t>> tmp2(n + m + 1,
                                        sequence<s_size_t>(n + m + 1, MAX_VAL));
      merge_horizontal(left, right, Matrix(tmp2, 0, 0), theta, n + m1 + 1,
                       n + m2 + 1, n);
      for (size_t x = 0; x < n + m + 1; x++) {
        for (size_t y = 0; y < n + m + 1; y++) {
          ret[x][y] = tmp2[x][y];
        }
      }
    } else if (m == 1) {
      auto up = Matrix(ret, 0, 0);
      auto down = Matrix(ret, n1 * 3, 0);
      solve_r(i, n1, j, m, r, c);
      solve_r(i + n1, n2, j, m, r + n1 * 3, c);
      sequence<sequence<s_size_t>> tmp2(n + m + 1,
                                        sequence<s_size_t>(n + m + 1, MAX_VAL));
      merge_vertical(up, down, Matrix(tmp2, 0, 0), theta, n1 + m + 1,
                     n2 + m + 1, m);
      for (size_t x = 0; x < n + m + 1; x++) {
        for (size_t y = 0; y < n + m + 1; y++) {
          ret[x][y] = tmp2[x][y];
        }
      }
    } else {
      size_t N = max(n1, m1) * 3;
      // printf("N: %zu\n", N);
      auto upper_left = Matrix(ret, 0, 0);
      auto bottom_left = Matrix(ret, N, 0);
      auto upper_right = Matrix(ret, 0, N);
      auto bottom_right = Matrix(ret, N, N);
      solve_r(i, n1, j, m1, r, c);
      solve_r(i + n1, n2, j, m1, r + N, c);
      // printf("here, i: %zu, n: %zu, j: %zu, m: %zu, r: %zu, c: %zu\n", i, n1,
      // j + m1, m2, r, c + N);
      solve_r(i, n1, j + m1, m2, r, c + N);
      solve_r(i + n1, n2, j + m1, m2, r + N, c + N);

      sequence<sequence<s_size_t>> tmp2(
          n + m1 + 1, sequence<s_size_t>(n + m1 + 1, MAX_VAL));
      sequence<sequence<s_size_t>> tmp3(
          n + m2 + 1, sequence<s_size_t>(n + m2 + 1, MAX_VAL));
      merge_vertical(upper_left, bottom_left, Matrix(tmp2, 0, 0), theta,
                     n1 + m1 + 1, n2 + m1 + 1, m1);
      // printf("i: %zu, n: %zu, j: %zu, m: %zu\n", i, n, j, m1);
      // for (size_t x = 0; x < n + m1 + 1; x++) {
      // for (size_t y = 0; y < n + m1 + 1; y++) {
      // printf("%10u%c", tmp2[x][y], " \n"[y == n + m1]);
      //}
      //}
      // puts("");

      merge_vertical(upper_right, bottom_right, Matrix(tmp3, 0, 0),
                     Matrix(theta, 0, N), n1 + m2 + 1, n2 + m2 + 1, m2);
      // printf("i: %zu, n: %zu, j: %zu, m: %zu\n", i, n, j + m1, m2);
      // for (size_t x = 0; x < n + m2 + 1; x++) {
      // for (size_t y = 0; y < n + m2 + 1; y++) {
      // printf("%10u%c", tmp3[x][y], " \n"[y == n + m2]);
      //}
      //}
      // puts("");
      merge_horizontal(Matrix(tmp2, 0, 0), Matrix(tmp3, 0, 0), ret, theta,
                       n + m1 + 1, n + m2 + 1, n);
    }
    // printf("i: %zu n: %zu, j: %zu m: %zu, r: %zu, c: %zu, use %zu by %zu\n",
    // i, n, j, m, r, c, n + m + 1, n + m + 1);
    // printf("i: %zu, n: %zu, j: %zu, m: %zu\n", i, n, j, m);
    // for (size_t x = 0; x < n + m + 1; x++) {
    // for (size_t y = 0; y < n + m + 1; y++) {
    // printf("%10u%c", ret[x][y], " \n"[y == n + m]);
    //}
    //}
    // puts("");
  }

  size_t solve() {
    size_t n = A.size(), m = B.size();
    size_t N = static_cast<size_t>(1) << log2_up(max(n, m));
    dist =
        sequence<sequence<s_size_t>>(N * 3, sequence<s_size_t>(N * 3, MAX_VAL));
    g_theta = sequence<sequence<s_size_t>>(
        N * 3, sequence<s_size_t>::uninitialized(N * 3));
    solve_r(0, n, 0, m, 0, 0);
    return dist[n][m];
  }
};

template class DAC_MM<uint32_t>;
