#ifndef SPARSE_TABLE_SEQUENTIAL_H_
#define SPARSE_TABLE_SEQUENTIAL_H_

template <typename Seq>
struct SparseTableSequential {
  unsigned int n;
  std::vector<unsigned int> log;
  std::vector<std::vector<typename Seq::value_type>> f;

  SparseTableSequential(const Seq& a) {
    n = a.size();
    log.resize(n + 1);
    log[1] = 0;
    for (unsigned int i = 2, x = 0; i <= n; i++) {
      if (i == (i & -i)) {
        log[i] = log[i - 1] + 1;
      } else {
        log[i] = log[i - 1];
      }
    }
    f.resize(n);
    for (unsigned int i = 0; i < n; i++) {
      f[i].resize(log[n] + 1);
      f[i][0] = a[i];
    }
    for (unsigned int j = 1; j <= log[n]; j++) {
      for (unsigned int i = 0; i + (1 << (j - 1)) < n; i++) {
        f[i][j] = std::min(f[i][j - 1], f[i + (1 << (j - 1))][j - 1]);
      }
    }
  }

  typename Seq::value_type Query(unsigned int l, unsigned int r) {
    unsigned int len = r - l + 1;
    unsigned int t = log[len];
    return std::min(f[l][t], f[r - (1 << t) + 1][t]);
  }
};

#endif  // SPARSE_TABLE_SEQUENTIAL_H_