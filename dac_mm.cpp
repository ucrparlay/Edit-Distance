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

  size_t get_pow2(size_t x) {
    size_t ret = 1;
    while ((ret << 1) < x) {
      ret <<= 1;
    }
    return ret;
  }

  sequence<sequence<s_size_t>> merge_horizontal(
      sequence<sequence<s_size_t>> &left, sequence<sequence<s_size_t>> &right,
      size_t k) {
    size_t n1 = left.size(), n2 = right.size();
    size_t n = n1 + n2 - k - 1;
    sequence<sequence<s_size_t>> ret(n, sequence<s_size_t>(n, MAX_VAL));
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

    sequence<sequence<s_size_t>> theta(n1, sequence<s_size_t>(n2, MAX_VAL));
    auto compute = [&](size_t i, size_t j, size_t l, size_t r) {
      if (i + k < j) {
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
      sequence<sequence<s_size_t>> &up, sequence<sequence<s_size_t>> &down,
      size_t k) {
    size_t n1 = up.size(), n2 = down.size();
    size_t n = n1 + n2 - k - 1;
    sequence<sequence<s_size_t>> ret(n, sequence<s_size_t>(n, MAX_VAL));
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

  sequence<sequence<s_size_t>> solve_r(size_t i, size_t n, size_t j, size_t m) {
    size_t n1 = n / 2, n2 = n - n1;
    size_t m1 = m / 2, m2 = m - m1;
    if (n == 1 && m == 1) {
      auto ret =
          sequence<sequence<s_size_t>>(3, sequence<s_size_t>(3, MAX_VAL));
      ret[0][0] = ret[2][2] = 0;
      ret[0][1] = ret[1][0] = ret[1][2] = ret[2][1] = 1;
      ret[1][1] = (A[i] != B[j]);
      return ret;
    } else if (n == 1) {
      sequence<sequence<s_size_t>> v1, v2;
      v1 = solve_r(i, n, j, m1);
      v2 = solve_r(i, n, j + m1, m2);
      auto ret = merge_horizontal(v1, v2, n);
      return ret;
    } else if (m == 1) {
      sequence<sequence<s_size_t>> v1, v2;
      v1 = solve_r(i, n1, j, m);
      v2 = solve_r(i + n1, n2, j, m);
      auto ret = merge_vertical(v1, v2, m);
      return ret;
    } else {
      sequence<sequence<s_size_t>> t1, t2, t3, t4;
      t1 = solve_r(i, n1, j, m1);
      t2 = solve_r(i + n1, n2, j, m1);
      t3 = solve_r(i, n1, j + m1, m2);
      t4 = solve_r(i + n1, n2, j + m1, m2);

      sequence<sequence<s_size_t>> v1, v2;
      v1 = merge_vertical(t1, t2, m1);
      v2 = merge_vertical(t3, t4, m2);

      auto ret = merge_horizontal(v1, v2, n);
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
