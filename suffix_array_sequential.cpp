#include "suffix_array_sequential.h"

#include <algorithm>
#include <array>
#include <vector>

namespace {

void Prepare(std::vector<int>& a) {
  std::vector<int> b = a;
  std::sort(b.begin(), b.end());
  for (int i = 0; i < a.size(); i++) {
    a[i] = lower_bound(b.begin(), b.end(), a[i]) - b.begin();
  }
}

}  // namespace

void SuffixArraySequential::Resort(std::vector<std::array<int, 2>>& key) {
  for (int i = 0; i < n; i++) key[i][0]++, key[i][1]++;
  std::vector<int> b(n + 1), sum(n + 2), at(n);
  for (int i = 0; i < n; i++) at[i] = b[key[i][1]]++;
  for (int i = 0; i <= n; i++) sum[i + 1] = sum[i] + b[i];
  for (int i = 0; i < n; i++) sa[sum[key[i][1]] + at[i]] = i;
  for (int i = 0; i <= n; i++) b[i] = 0;
  for (int i = 0; i < n; i++) at[sa[i]] = b[key[sa[i]][0]]++;
  for (int i = 0; i <= n; i++) sum[i + 1] = sum[i] + b[i];
  for (int i = 0; i < n; i++) sa[sum[key[i][0]] + at[i]] = i;
  for (int i = 0, r = 0; i < n; i++) {
    if (i > 0 && key[sa[i]] > key[sa[i - 1]]) r++;
    rank[sa[i]] = r;
  }
}

void SuffixArraySequential::Build(const std::vector<int>& a_) {
  n = a_.size();
  a.resize(n);
  for (int i = 0; i < n; i++) a[i] = a_[i];
  Prepare(a);
  a.resize(n * 2, -1);
  sa.resize(n);
  rank.resize(n);
  std::vector<std::array<int, 2>> key(n);
  for (int i = 0; i < n; i++) key[i] = {a[i], -1};
  Resort(key);
  for (int len = 1; len <= n; len *= 2) {
    for (int i = 0; i < n; i++) {
      int r1 = rank[i];
      int r2 = i + len < n ? rank[i + len] : -1;
      key[i] = {r1, r2};
    }
    Resort(key);
  }
  height.resize(n);
  for (int i = 0, p = 0; i < n; i++) {
    if (rank[i] == 0) continue;
    int j = sa[rank[i] - 1];
    while (a[i + p] == a[j + p]) p++;
    height[rank[i]] = p;
    if (p) p--;
  }
}
