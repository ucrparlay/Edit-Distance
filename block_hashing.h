#ifndef BLOCK_HASHING_H
#define BLOCK_HASHING_H
#include "parlay/primitives.h"
#include "parlay/sequence.h"

using namespace parlay;
using namespace std;
using Hash_Type = uint32_t;

static constexpr size_t PRIME = 479;

uint32_t qpow2(size_t n) {
  uint32_t ret = 1;
  uint32_t a = PRIME;
  while (n) {
    if (n & 1) {
      ret = ret * a;
    }
    a = a * a;
    n >>= 1;
  }
  return ret;
}

template <typename Seq, typename s_size_t = uint32_t>
class Block_Hashing {
  static constexpr int LOG2_BLOCK = 6;
  static constexpr size_t BLOCK_SIZE = 1 << LOG2_BLOCK;
  static constexpr size_t BLOCK_MASK = BLOCK_SIZE - 1;
  static constexpr size_t GRANULARITY = 1024;
  const Seq& A;
  const Seq& B;
  sequence<sequence<s_size_t>> blocking_A;
  sequence<sequence<s_size_t>> blocking_B;
  sequence<pair<Hash_Type, Hash_Type>> pre_sub_A;
  sequence<pair<Hash_Type, Hash_Type>> pre_sub_B;

 public:
  Block_Hashing(const Seq& _A, const Seq& _B) : A(_A), B(_B) {
    blocking_A = build_blocking(A);
    blocking_B = build_blocking(B);
    pre_sub_A = build_pre_sub(A);
    pre_sub_B = build_pre_sub(B);
  }
  sequence<sequence<s_size_t>> build_blocking(const Seq& S) {
    size_t n = (S.size() + BLOCK_SIZE - 1) / BLOCK_SIZE;
    size_t k = max(size_t{1}, (size_t)ceil(log2(n)));
    auto blocking =
        sequence<sequence<s_size_t>>(k, sequence<Hash_Type>::uninitialized(n));
    parallel_for(0, n, [&](size_t i) {
      size_t s = i << LOG2_BLOCK;
      size_t e = min(S.size(), (i + 1) << LOG2_BLOCK);
      Hash_Type p = 1;
      blocking[0][i] = 0;
      for (size_t j = 0; j < e - s; j++) {
        blocking[0][i] += S[e - 1 - j] * p;
        p *= PRIME;
      }
    });
    for (size_t i = 1; i < k; i++) {
      parallel_for(
          0, n - (1 << i) + 1,
          [&](size_t j) {
            blocking[i][j] =
                blocking[i - 1][j] * qpow2((1 << (i - 1)) << LOG2_BLOCK) +
                blocking[i - 1][j + (1 << (i - 1))];
          },
          GRANULARITY);
    }
    // printf("k: %zu, n: %zu\n", k, n);
    // for (size_t i = 0; i < k; i++) {
    // printf("blocking[%zu]: ", i);
    // for (size_t j = 0; j < n - (1 << i) + 1; j++) {
    // printf("%u ", blocking[i][j]);
    //}
    // puts("");
    //}
    return blocking;
  }
  sequence<pair<Hash_Type, Hash_Type>> build_pre_sub(const Seq& S) {
    size_t n = S.size();
    size_t num_blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    auto pre_sub = sequence<pair<Hash_Type, Hash_Type>>::uninitialized(n);
    parallel_for(0, num_blocks, [&](size_t i) {
      size_t s = i << LOG2_BLOCK;
      size_t e = min(n, (i + 1) << LOG2_BLOCK);
      pre_sub[s].first = S[s];
      for (size_t j = 0; j < e - s - 1; j++) {
        pre_sub[s + j + 1].first = pre_sub[s + j].first * PRIME + S[s + j + 1];
      }
      pre_sub[e - 1].second = S[e - 1];
      Hash_Type p = PRIME;
      for (size_t j = 0; j < e - s - 1; j++) {
        pre_sub[e - j - 2].second =
            pre_sub[e - j - 1].second + p * S[e - j - 2];
        p *= PRIME;
      }
    });
    // for (size_t i = 0; i < n; i++) {
    // printf("pre_sub[%zu]: (%u,%u)\n", i, pre_sub[i].first,
    // pre_sub[i].second);
    //}
    return pre_sub;
  }
  Hash_Type get_value(const sequence<sequence<Hash_Type>>& blocking,
                      const sequence<pair<Hash_Type, Hash_Type>>& pre_sub,
                      size_t pos, size_t len, size_t log2_len) {
    if (!(pos & BLOCK_MASK)) {
      size_t block_id = pos >> LOG2_BLOCK;
      return blocking[log2_len][block_id] * qpow2(BLOCK_SIZE) +
             blocking[0][block_id + (1 << log2_len)];
    } else {
      size_t tail = pos & BLOCK_MASK;
      size_t head = BLOCK_SIZE - tail;
      // size_t head = BLOCK_SIZE - pos % BLOCK_SIZE;
      // size_t tail = (pos + len) % BLOCK_SIZE;
      assert(head != 0 && tail != 0 && head + tail == BLOCK_SIZE);

      size_t block_id = (pos + head) >> LOG2_BLOCK;
      // Combine head
      size_t ret = pre_sub[pos].second;
      // Combine blocks;
      ret = ret * qpow2((1 << log2_len) << LOG2_BLOCK) +
            blocking[log2_len][block_id];
      // Combine tail
      ret = ret * qpow2(tail) + pre_sub[pos + len - 1].first;
      return ret;
    }
  }
  int query_lcp(size_t x, size_t y) {
    auto xx = x, yy = y;
    int log2_len = 0;
    size_t block_len = (1 << log2_len) * BLOCK_SIZE;
    size_t len = block_len + BLOCK_SIZE;
    size_t ret = 0;
    auto if_same = [&]() {
      if (x + len > A.size() || y + len > B.size()) {
        return false;
      }
      return get_value(blocking_A, pre_sub_A, x, len, log2_len) ==
             get_value(blocking_B, pre_sub_B, y, len, log2_len);
    };
    while (true) {
      if (if_same()) {
        ret = len;
        log2_len++;
        block_len <<= 1;
        len = block_len + BLOCK_SIZE;
      } else {
        break;
      }
    }
    x += ret;
    y += ret;
    auto tmp = log2_len;
    if (log2_len) {
      do {
        log2_len--;
        block_len >>= 1;
        len = block_len + BLOCK_SIZE;
        if (if_same()) {
          ret += len;
          x += len;
          y += len;
        }
      } while (log2_len);
    }
    int cnt = 0;
    while (x < A.size() && y < B.size() && A[x] == B[y]) {
      ret++;
      x++;
      y++;
      cnt++;
    }
    // if (cnt >= 100) {
    // printf("tmp: %zu\n", tmp);
    // printf("x: %zu, y: %zu, ret: %d\n", xx, yy, ret);
    // printf("block_size: %zu\n", BLOCK_SIZE);
    // log2_len = 0;
    // block_len = (1 << log2_len) * BLOCK_SIZE;
    // len = block_len + BLOCK_SIZE;
    // printf("len: %zu, v1: %u, v2: %u\n", len,
    // get_value(blocking_A, pre_sub_A, xx, len, log2_len),
    // get_value(blocking_B, pre_sub_B, yy, len, log2_len));
    // printf("log2_len: %zu\n", log2_len);
    // if (xx % BLOCK_SIZE == 0) {
    // printf("v1: %u * %u + %u\n", blocking_A[log2_len][xx],
    // qpow2(BLOCK_SIZE), blocking_A[0][xx + (1 << log2_len)]);
    //} else {
    // size_t tail = xx & BLOCK_MASK;
    // size_t head = BLOCK_SIZE - tail;
    // size_t block_id = (xx + head) >> LOG2_BLOCK;
    // printf("v1: %u\n", pre_sub_A[xx].second);
    // printf("v1 = v1 * %u + %u\n", qpow2((1 << log2_len) << LOG2_BLOCK),
    // blocking_A[log2_len][block_id]);
    // printf("v1 = v1 * %u + %u\n", qpow2(tail),
    // pre_sub_A[xx + len - 1].first);
    //}
    // if (yy % BLOCK_SIZE == 0) {
    // printf("v2: %u * %u + %u\n", blocking_B[log2_len][yy],
    // qpow2(BLOCK_SIZE), blocking_B[0][yy + (1 << log2_len)]);
    //} else {
    // size_t tail = yy & BLOCK_MASK;
    // size_t head = BLOCK_SIZE - tail;
    // size_t block_id = (yy + head) >> LOG2_BLOCK;
    // printf("v2: %u\n", pre_sub_B[yy].second);
    // printf("v2 = v2 * %u + %u\n", qpow2((1 << log2_len) << LOG2_BLOCK),
    // blocking_B[log2_len][block_id]);
    // printf("v2 = v2 * %u + %u\n", qpow2(tail),
    // pre_sub_B[yy + len - 1].first);
    //}
    // log2_len = 0;
    // for (size_t i = 0; i < 10 && xx + i < A.size(); i++) {
    // printf("A[%zu]: %u\n", xx + i, A[xx + i]);
    //}
    // printf("\n");
    // for (size_t i = 0; i < 10 && yy + i < B.size(); i++) {
    // printf("B[%zu]: %u\n", yy + i, B[yy + i]);
    //}
    //}
    assert(cnt < 2 * BLOCK_SIZE);
    // assert(xx + ret <= A.size() && yy + ret <= B.size());
    // for (size_t i = 0; i < ret; i++) {
    // if (A[xx + i] != B[yy + i]) {
    // printf("tmp: %zu\n", tmp);
    // printf("x: %zu, y: %zu, ret: %d\n", xx, yy, ret);
    // printf("block_size: %zu\n", BLOCK_SIZE);
    // log2_len = 0;
    // block_len = (1 << log2_len) * BLOCK_SIZE;
    // len = block_len + BLOCK_SIZE;
    // printf("len: %zu, v1: %u, v2: %u\n", len,
    // get_value(blocking_A, pre_sub_A, xx, len, log2_len),
    // get_value(blocking_B, pre_sub_B, yy, len, log2_len));
    // printf("log2_len: %zu\n", log2_len);
    // if (xx % BLOCK_SIZE == 0) {
    // printf("v1: %u * %u + %u\n", blocking_A[log2_len][xx],
    // qpow2(BLOCK_SIZE), blocking_A[0][xx + (1 << log2_len)]);
    //} else {
    // size_t tail = xx & BLOCK_MASK;
    // size_t head = BLOCK_SIZE - tail;
    // size_t block_id = (xx + head) >> LOG2_BLOCK;
    // printf("v1: %u\n", pre_sub_A[xx].second);
    // printf("v1 = v1 * %u + %u\n", qpow2((1 << log2_len) << LOG2_BLOCK),
    // blocking_A[log2_len][block_id]);
    // printf("v1 = v1 * %u + %u\n", qpow2(tail),
    // pre_sub_A[xx + len - 1].first);
    //}
    // if (yy % BLOCK_SIZE == 0) {
    // printf("v2: %u * %u + %u\n", blocking_B[log2_len][yy],
    // qpow2(BLOCK_SIZE), blocking_B[0][yy + (1 << log2_len)]);
    //} else {
    // size_t tail = yy & BLOCK_MASK;
    // size_t head = BLOCK_SIZE - tail;
    // size_t block_id = (yy + head) >> LOG2_BLOCK;
    // printf("v2: %u\n", pre_sub_B[yy].second);
    // printf("v2 = v2 * %u + %u\n", qpow2((1 << log2_len) << LOG2_BLOCK),
    // blocking_B[log2_len][block_id]);
    // printf("v2 = v2 * %u + %u\n", qpow2(tail),
    // pre_sub_B[yy + len - 1].first);
    //}
    // log2_len = 0;
    // for (size_t i = 0; i < 10 && xx + i < A.size(); i++) {
    // printf("A[%zu]: %u\n", xx + i, A[xx + i]);
    //}
    // printf("\n");
    // for (size_t i = 0; i < 10 && yy + i < B.size(); i++) {
    // printf("B[%zu]: %u\n", yy + i, B[yy + i]);
    //}
    //}
    // assert(A[xx + i] == B[yy + i]);
    //}
    return ret;
  }
};

#endif
