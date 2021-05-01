//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  ppval.cpp
//
//  Code generation for function 'ppval'
//


// Include files
#include "ppval.h"
#include "rt_nonfinite.h"
#include "test_controller.h"
#include <cstring>

// Function Definitions
void ppval(const double pp_breaks[100], const double pp_coefs[3564], double x,
           double v[9])
{
  if (rtIsNaN(x)) {
    for (int mid_i = 0; mid_i < 9; mid_i++) {
      v[mid_i] = x;
    }
  } else {
    int low_i;
    int mid_i;
    int low_ip1;
    int high_i;
    int icp;
    double xloc;
    low_i = 1;
    low_ip1 = 2;
    high_i = 100;
    while (high_i > low_ip1) {
      mid_i = (low_i + high_i) >> 1;
      if (x >= pp_breaks[mid_i - 1]) {
        low_i = mid_i;
        low_ip1 = mid_i + 1;
      } else {
        high_i = mid_i;
      }
    }

    icp = (low_i - 1) * 9;
    xloc = x - pp_breaks[low_i - 1];
    std::memcpy(&v[0], &pp_coefs[icp], 9U * sizeof(double));
    for (low_ip1 = 0; low_ip1 < 3; low_ip1++) {
      high_i = icp + (low_ip1 + 1) * 891;
      for (mid_i = 0; mid_i < 9; mid_i++) {
        v[mid_i] = xloc * v[mid_i] + pp_coefs[high_i + mid_i];
      }
    }
  }
}

// End of code generation (ppval.cpp)
