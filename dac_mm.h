#include <parlay/delayed_sequence.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

using namespace parlay;

template <typename T, typename s_size_t = uint32_t>
class DAC_MM {
 protected:
  const sequence<T> &A;
  const sequence<T> &B;

  size_t get_pow2(size_t);
  size_t sqr(size_t);
  size_t get_size(size_t, size_t);

  sequence<sequence<uint32_t>> merge_horizontal(
      sequence<sequence<uint32_t>> &left, sequence<sequence<uint32_t>> &right,
      size_t k);

  sequence<sequence<uint32_t>> merge_vertical(
      sequence<sequence<uint32_t>> &up, sequence<sequence<uint32_t>> &down,
      size_t k);

  void print_matrix(size_t i, size_t n, size_t j, size_t m,
                    sequence<sequence<uint32_t>> &ret);

  sequence<sequence<uint32_t>> solve_r(size_t i, size_t n, size_t j, size_t m);

 public:
  DAC_MM(const sequence<T> &_A, const sequence<T> &_B);
  size_t solve();
};
