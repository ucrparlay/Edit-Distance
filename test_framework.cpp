#include <parlay/delayed_sequence.h>
#include <parlay/internal/get_time.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/utilities.h>

#include <fstream>

#include "dac_mm.h"
#include "dac_mm_k.h"
#include "edit_distance_block_hashing.h"
#include "edit_distance_dp.h"
#include "edit_distance_hashing.h"
#include "edit_distance_parallel.h"
#include "minimum_edit_distance.h"

constexpr size_t NUM_TESTS = 6;
size_t num_rounds = 3;

template <typename T>
auto generate_strings(size_t n, size_t k, size_t alpha, size_t seed = 0) {
  printf("Generating test case... (n: %zu, k: %zu, alpha: %zu)\n", n, k, alpha);
  parlay::sequence<T> A(n), B(n);
  parlay::parallel_for(
      0, n, [&](size_t i) { A[i] = B[i] = parlay::hash32(i + seed) % alpha; });

  // substitions, insertions, and deletions and roughly equally distributed
  size_t _k = k / 3;

  // substitutions
  parlay::parallel_for(0, _k, [&](size_t i) {
    size_t idx = parlay::hash32(i + seed) % n;
    B[idx] = parlay::hash32(i + n + seed) % alpha;
  });

  // insertions and deletions
  auto pred1 = parlay::delayed_seq<bool>(
      n, [&](size_t i) { return parlay::hash32_2(i + seed) % n >= _k; });
  auto pred2 = parlay::delayed_seq<bool>(
      n, [&](size_t i) { return parlay::hash32_2(i + n + seed) % n >= _k; });
  A = pack(A, pred1);
  B = pack(B, pred2);
  return std::make_tuple(A, B);
}

std::string test_name(int id) {
  switch (id) {
    case 0:
      return "dp";
      break;
    case 1:
      return "string_hashing_bfs";
      break;
    case 2:
      return "parallel bfs";
      break;
    case 3:
      return "parlay::edit_distance";
      break;
    case 4:
      return "dac_mm";
      break;
    case 5:
      return "dac_mm_k";
      break;
    case 6:
      return "block_hashing";
      break;
    default:
      assert(0);
  }
}

template <typename T>
double test(const parlay::sequence<T> &A, const parlay::sequence<T> &B,
            int id) {
  std::cout << "\nTest name: " << test_name(id) << std::endl;
  double total_time = 0;
  for (size_t i = 0; i <= num_rounds; i++) {
    parlay::internal::timer t;
    size_t num_edits;
    switch (id) {
      case 0:
        num_edits = EditDistanceDP<T>().Solve(A, B);
        break;
      case 1:
        num_edits = EditDistanceHashParallel(A, B);
        break;
      case 2:
        num_edits = EditDistanceParallel().Solve(A, B);
        break;
      case 3:
        num_edits = minimum_edit_distance(A, B);
        break;
      case 4:
        num_edits = DAC_MM<sequence<uint32_t>>(A, B).solve();
        break;
      case 5:
        num_edits = DAC_MM_K<sequence<uint32_t>>(A, B).solve();
        break;
      case 6:
        num_edits = EditDistanceBlockHashParallel(A, B);
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
  std::ofstream ofs("edit_distance.tsv", std::ios_base::app);
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
        "rounds: number of rounds");
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
  // for (size_t i = 1; i <= 500; i++) {
  // for (size_t j = 3 * i; j >= 1; j -= 3) {
  // for (size_t seed = 0; seed < 1; seed++) {
  // printf("i: %zu, j: %zu\n", i, j);
  // parlay::sequence<Type> A, B;
  // std::tie(A, B) = generate_strings<Type>(i, j, 3 * i, seed);
  ////printf("A.size(): %zu, B.size(): %zu\n", A.size(), B.size());
  ////printf("A: ");
  ////for (size_t k = 0; k < A.size(); k++) {
  ////printf("%u ", A[k]);
  ////}
  ////puts("");
  ////printf("B: ");
  ////for (size_t k = 0; k < B.size(); k++) {
  ////printf("%u ", B[k]);
  ////}
  ////puts("");
  // size_t v1 = EditDistanceDP<Type>().Solve(A, B);
  // size_t v2 = DAC_MM_K<sequence<Type>>(A, B).solve();
  // if (v1 != v2) {
  // printf("v1: %zu, v2: %zu\n", v1, v2);
  // printf("wrong answer\n");
  // if (A.size() < 20) {
  // return 0;
  //} else {
  // getchar();
  //}
  //}
  //}
  //}
  //}
  parlay::sequence<Type> A, B;
  std::tie(A, B) = generate_strings<Type>(n, k, alpha);
  run_all(A, B, 5);

  /*
    for real datasets
  */
  // std::string str_A;
  // std::string str_B;
  // parse_text_file_with_blank("./data_prep/data/1.txt", A);
  // parse_text_file_with_blank("./data_prep/data/2.txt", B);
  // printf("size A: %d\n", A.size());
  // printf("size B: %d\n", B.size());
  // run_all(A, B, 0);
  // run_all(A, B, 1);
  // run_all(A, B, 2);
  // run_all(A, B, 5);
  // run_all(A, B, 3);

  return 0;
}
