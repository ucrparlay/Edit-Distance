//
//  utils.h
//  ed
//
//  Created by Alan on 11/27/22.
//

#ifndef utils_h
#define utils_h

#include <malloc.h>
#include <math.h>
#include <stdint.h>
#include <sys/resource.h>
#include <unistd.h>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "parlay/internal/get_time.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

constexpr uint32_t PRIME_BASE = 479;
constexpr uint32_t PRIME = 479;

// const int PRIME_BASE = 631;

#define FASTLOG2(X) \
  ((unsigned)(8 * sizeof(unsigned long long) - __builtin_clzll((X)) - 1))

// /**
//  * Inplace block scan DIRECTLY without filter
//  */
// template <typename Seq>
// void s_inplace_scan_inclusive(Seq& A, size_t n) {
//   auto block_size = std::max((size_t)(8000), (size_t)std::sqrt(r - l));
//   if (r - l <= 100000) {
//     for (size_t i = 1; i < n; i++) {
//       A[i] += A[i - 1] * PRIME;
//     }
//   } else {
//     size_t num_blocks = (r - l - 1) / block_size + 1;
//     parlay::parallel_for(0, num_blocks, [&](size_t i) {
//       // aux[i] = inplace_seq_scan_exclusive_direct(
//       //     A, l + i * block_size, std::min(l + i * block_size + block_size,
//       //     r));
//       for (size_t j = i * block_size + 1; j < std::min((i + 1) * block_size,
//       n);
//            j++) {
//         A[i] += A[i - 1] * PRIME;
//       }
//     });

//     for (size_t k = 1; k < num_blocks; k++) {
//       A[(k + 1) * block_size - 1] +=
//           A[k * block_size - 1] * quick_power(block_size);
//     }

//     parlay::parallel_for(1, num_blocks, [&](size_t i) {
//       for (size_t j = i * block_size;
//            j < std::min(i * block_size + block_size, r); j++) {
//         A[j] += A[block_size * (j / block_size) - 1];
//       }
//     });
//   }
// }

template <typename T>
std::vector<T> get_unique_elems(const std::vector<T>& v) {
  std::vector<T> unique_elems;
  std::unordered_set<T> seen_elems;
  for (const auto& elem : v) {
    if (seen_elems.count(elem) == 0) {
      seen_elems.insert(elem);
      unique_elems.push_back(elem);
    }
  }
  return unique_elems;
}

// Function for mapping chars to bit-wise indices
template <typename T>
parlay::sequence<int> map_to_indices(const std::vector<T>& vec) {
  parlay::sequence<int> indices(vec.size());
  std::vector<T> unique_elems = get_unique_elems(vec);
  std::unordered_map<T, int> elem_to_index;
  for (int i = 0; i < unique_elems.size(); i++) {
    elem_to_index[unique_elems[i]] = i;
  }
  for (int i = 0; i < vec.size(); i++) {
    indices[i] = elem_to_index[vec[i]];
  }
  return indices;
}

// Function to parse a file with single-line string to
// parlay sequence.
// Usage: `parse_string_line(file_path, destination_variables)`
template <typename T>
void parse_string_line(const std::string file_name,
                       parlay::sequence<T>& parse_var) {
  std::ifstream file(file_name.c_str());
  if (!file) {
    printf("Cannot open file \n");
    return;
  }
  std::string str;
  std::getline(file, str);
  for (auto iter : str) {
    parse_var.push_back((uint32_t)iter);
  }
}

// Function to parse a file with a text file to
// parlay::sequence variable, blank lines and
// white spaces will be omitted
// Usage: `parse_text_file_with_blank(file_path, destination_variables)`
template <typename T>
void parse_text_file_with_blank(const std::string filename,
                                parlay::sequence<T>& parser_var) {
  char ch;
  std::fstream fin(filename, std::fstream::in);
  while (fin >> std::noskipws >> ch) {
    unsigned char valid_ch = (unsigned char)(ch);
    parser_var.push_back((T)valid_ch);
  }
}

class Timer {
 public:
  Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }

  double elapsed() {
    clock_gettime(CLOCK_REALTIME, &end_);
    return end_.tv_sec - beg_.tv_sec +
           (end_.tv_nsec - beg_.tv_nsec) / 1000000000.;
  }

  void reset() { clock_gettime(CLOCK_REALTIME, &beg_); }

 private:
  timespec beg_, end_;
};

// // /**
// //  * Returns the peak (maximum so far) resident set size (physical
// //  * memory use) measured in bytes, or zero if the value cannot be
// //  * determined on this OS.
// //  */
// size_t getPeakRSS() {
//   /* BSD, Linux, and OSX -------------------------------------- */
//   struct rusage rusage;
//   getrusage(RUSAGE_SELF, &rusage);
//   return (size_t)(rusage.ru_maxrss * 1024L);
// }

#endif /* utils_h */
