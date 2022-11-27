#include <parlay/delayed_sequence.h>
#include <parlay/internal/get_time.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/utilities.h>

#include <fstream>

#include "minimum_edit_distance.h"
constexpr size_t NUM_TESTS = 1;
size_t num_rounds = 10;

template <typename T>
auto generate_strings(size_t n, size_t k, size_t alpha) {
  printf("Generating test case... (n: %zu, k: %zu, alpha: %zu)\n", n, k, alpha);
  parlay::sequence<T> A(n), B(n);
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

std::string test_name(int id) {
  switch (id) {
    case 0:
      return "parlay::edit_distance";
    default:
      assert(0);
  }
}

template <typename T>
double test(const parlay::sequence<T> &A, const parlay::sequence<T> &B,
            int id) {
  std::cout << "Test name: " << test_name(id) << std::endl;
  double total_time = 0;
  for (size_t i = 0; i <= num_rounds; i++) {
    parlay::internal::timer t;
    size_t num_edits;
    switch (id) {
      case 0:
        num_edits = minimum_edit_distance(A, B);
        break;
      default:
        assert(0);
    }
    t.stop();
    if (i == 0) {
      printf("#edits: %zu\n", num_edits);
      printf("Warmup round: %f\n", t.total_time());
    } else {
      printf("Round %zu: %f\n", i, t.total_time());
      total_time += t.total_time();
    }
  }
  double average_time = total_time / num_rounds;
  printf("Average time: %f\n", total_time / num_rounds);
  return average_time;
}

template <typename T>
void run_all(const parlay::sequence<T> &A, const parlay::sequence<T> &B,
             int id = -1) {
  std::vector<double> times;
  if (id == -1) {
    for (size_t i = 0; i < NUM_TESTS; i++) {
      times.push_back(test(A, B, i));
    }
  } else {
    times.push_back(test(A, B, id));
  }
  std::ofstream ofs("edit_distance.tsv");
  for (auto t : times) {
    ofs << t << '\t';
  }
  ofs << '\n';
  ofs.close();
}

int main(int argc, char *argv[]) {
  size_t n = 1000000;
  size_t k = 1000;
  size_t alpha = n;
  if (argc == 1) {
    printf(
        "Usage: ./edit_distance <n> <k> <alpha> <rounds>\n"
        "n: length of strings\n"
        "k: estimated number of edits\n"
        "alpha: alphabet size\n"
        "rounsd: number of rounds");
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
  using Type = uint32_t;
  parlay::sequence<Type> A, B;
  std::tie(A, B) = generate_strings<Type>(n, k, alpha);
  run_all(A, B);

  return 0;
}
