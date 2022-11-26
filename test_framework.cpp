#include <parlay/delayed_sequence.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/utilities.h>

using Type = uint32_t;

auto generate_strings(size_t n, size_t k, size_t alpha) {
  parlay::sequence<Type> A(n), B(n);
  parlay::parallel_for(
      0, n, [&](size_t i) { A[i] = B[i] = parlay::hash32(i) % alpha; });

  // substitutions
  parlay::parallel_for(0, k, [&](size_t i) {
    size_t idx = parlay::hash32(i) % n;
    B[idx] = parlay::hash32(i + n) % alpha;
  });

  // insertions and deletions
  auto pred1 = parlay::delayed_seq<bool>(
      n, [&](size_t i) { return parlay::hash32_2(i) % n >= k; });
  auto pred2 = parlay::delayed_seq<bool>(
      n, [&](size_t i) { return parlay::hash32_2(i + n) % n >= k; });
  A = pack(A, pred1);
  B = pack(B, pred2);
  return std::make_tuple(A, B);
}

int main(int argc, char* argv[]) {
  size_t n = 1000000;
  size_t k = 1000;
  size_t alpha = n;
  if (argc == 1) {
    printf(
        "Usage: ./edit_distance <n> <k> <alpha>\n"
        "n: length of strings\n"
        "k: estimated number of edits\n"
        "alpha: alphabet size\n");
    exit(0);
  }
  if (argc >= 2) {
    n = atoi(argv[1]);
  }
  if (argc >= 3) {
    k = atoi(argv[2]);
  }
  if (argc >= 4) {
    alpha = atoi(argv[3]);
  }
  parlay::sequence<Type> A, B;
  std::tie(A, B) = generate_strings(n, k, alpha);
  // size_t edits = edit_distance(A, B);
  // printf("#edits: %zu\n", edits);
  return 0;
}
