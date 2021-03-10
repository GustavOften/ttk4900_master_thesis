//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  mod.cpp
//
//  Code generation for function 'mod'
//


// Include files
#include "mod.h"
#include "rt_nonfinite.h"
#include "test_controller.h"
#include <math.h>

// Function Definitions
double b_mod(double x)
{
  double r;
  if (rtIsNaN(x) || rtIsInf(x)) {
    r = rtNaN;
  } else if (x == 0.0) {
    r = 0.0;
  } else {
    bool rEQ0;
    r = fmod(x, 3.1415926535897931);
    rEQ0 = (r == 0.0);
    if (!rEQ0) {
      double q;
      q = fabs(x / 3.1415926535897931);
      rEQ0 = !(fabs(q - floor(q + 0.5)) > 2.2204460492503131E-16 * q);
    }

    if (rEQ0) {
      r = 0.0;
    } else {
      if (x < 0.0) {
        r += 3.1415926535897931;
      }
    }
  }

  return r;
}

// End of code generation (mod.cpp)
