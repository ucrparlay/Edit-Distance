//
//  utils.h
//  ed
//
//  Created by Alan on 11/27/22.
//

#ifndef utils_h
#define utils_h

#include <math.h>
#include <stdint.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "parlay/internal/get_time.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

const int PRIME_BASE = 379;
// const int PRIME_BASE = 631;

#define FASTLOG2(X) \
  ((unsigned)(8 * sizeof(unsigned long long) - __builtin_clzll((X)) - 1))

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
    parser_var.push_back((uint32_t)(ch));
  }
}

#endif /* utils_h */
