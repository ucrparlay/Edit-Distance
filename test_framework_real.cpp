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

constexpr size_t NUM_TESTS = 4;
size_t num_rounds = 3;

class InputParser {
 public:
  InputParser(int &argc, char **argv) {
    for (int i = 1; i < argc; ++i) this->tokens.push_back(std::string(argv[i]));
  }

  const std::string &getCmdOption(const std::string &option) const {
    std::vector<std::string>::const_iterator itr;
    itr = std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
      return *itr;
    }
    static const std::string empty_string("");
    return empty_string;
  }

  bool cmdOptionExists(const std::string &option) const {
    return std::find(this->tokens.begin(), this->tokens.end(), option) !=
           this->tokens.end();
  }

 private:
  std::vector<std::string> tokens;
};

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
      return "BFS-Hash";
      break;
    case 1:
      return "BFS-B-Hash";
      break;
    case 2:
      return "BFS-SA";
      break;
    case 3:
      return "DaC-MM-K";
      break;
    case 4:
      return "DaC-MM";
      break;
    case 5:
      return "DP";
      break;
    case 6:
      return "ParlayLib";
      break;
    default:
      assert(0);
  }
}

template <typename T>
double test(const parlay::sequence<T> &A, const parlay::sequence<T> &B,
            int id) {
  // std::cout << "\nTest name: " << test_name(id) << std::endl;
  double total_time = 0;
  double building_time_total = 0;
  for (size_t i = 0; i <= num_rounds; i++) {
    parlay::internal::timer t;
    double b_time;
    size_t num_edits;
    switch (id) {
      case 0:
        num_edits = EditDistanceHashParallel(A, B, &b_time);
        break;
      case 1:
        num_edits = EditDistanceBlockHashParallel(A, B, &b_time);
        break;
      case 2:
        num_edits = EditDistanceParallel().Solve(A, B, &b_time);
        break;
      case 3:
        num_edits = DAC_MM_K<sequence<uint32_t>>(A, B).solve();
        break;
      case 4:
        num_edits = DAC_MM<sequence<uint32_t>>(A, B).solve();
        break;
      case 5:
        num_edits = EditDistanceDP<T>().Solve(A, B);
        break;
      case 6:
        num_edits = minimum_edit_distance(A, B);
        break;
      default:
        assert(0);
    }
    t.stop();
    if (i == 0) {
      std::cout << test_name(id) << std::endl;
      printf("#edits: %zu\n", num_edits);
      printf("Warmup round: %f, %f\n", b_time, t.total_time());
    } else {
      printf("Round %zu: %f\n", i, t.total_time());
      total_time += t.total_time();
      building_time_total += b_time;
    }
  }
  double average_time = total_time / num_rounds;
  printf("avg building time: %f, ", building_time_total / num_rounds);
  printf("avg time: %f\n", total_time / num_rounds);

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
  int id = -1;
  std::string path_1;
  std::string path_2;

  InputParser input(argc, argv);
  id = atoi(argv[1]);

  const std::string &filename1 = input.getCmdOption("-f1");
  if (!filename1.empty()) {
    path_1 = filename1;
  }

  const std::string &filename2 = input.getCmdOption("-f2");
  if (!filename2.empty()) {
    path_2 = filename2;
  }
  if (argc == 1) {
    printf(
        "Usage: ./edit_distance -i <id> -n <n> -k <k> -a <alpha> -r <rounds> "
        "-f1 <file_path1> "
        " -f2 <file_path2>\n"
        "id: id of the algorithm\n"
        "n: length of strings\n"
        "k: estimated number of edits\n"
        "alpha: alphabet size\n"
        "rounds: number of rounds\n"
        "file_path1: text file 1\n"
        "file_path2: text file 2");
    exit(0);
  }

  /*
    for real datasets
  */
  using Type = uint32_t;
  parlay::sequence<Type> A, B;
  parse_text_file_with_blank(path_1, A);
  parse_text_file_with_blank(path_2, B);
  printf("size A: %d\n", A.size());
  printf("size B: %d\n", B.size());
  // run_all(A, B, 6);
  run_all(A, B, id);

  return 0;
}
