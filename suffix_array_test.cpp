#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "parlay/parallel.h"
#include "suffix_array_parallel.h"
#include "suffix_array_sequential.h"

using namespace std;

int main() {
  string s = "abcaa";
  auto sa1 = suffix_array(s);
  SuffixArraySequential sa_seq;
  sa_seq.Build(s);
  auto sa2 = sa_seq.sa;
  int n = s.length();
  assert(sa1.size() == sa2.size());
  for (int i = 0; i < n; i++) {
    assert(sa1[i] == sa2[i]);
  }
  cout << "Pass: suffix_array_test" << endl;
  return 0;
}
