#include <fstream>
#include <iostream>
#include <string>

#include "edit_distance_sequential.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: ./edit_distance_test <filename>\n";
    return 0;
  }
  EditDistanceSequential* ed = new EditDistanceSequential();
  std::string s1, s2;
  std::ifstream fin(argv[1]);
  fin >> s1 >> s2;
  std::cout << s1 << '\n' << s2 << std::endl;
  std::vector<int> a(s1.length()), b(s2.length());
  for (unsigned int i = 0; i < s1.length(); i++) a[i] = s1[i];
  for (unsigned int i = 0; i < s2.length(); i++) b[i] = s2[i];
  std::cout << ed->Solve(a, b) << std::endl;
}