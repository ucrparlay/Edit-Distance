#ifndef EDIT_DISTANCE_PARALLAL_H_
#define EDIT_DISTANCE_PARALLAL_H_

#include "parlay/sequence.h"

size_t EditDistanceSA(const parlay::sequence<uint32_t>& a,
                      const parlay::sequence<uint32_t>& b, double* building_tm,
                      bool use_DC3 = false);

#endif  // EDIT_DISTANCE_PARALLAL_H_
