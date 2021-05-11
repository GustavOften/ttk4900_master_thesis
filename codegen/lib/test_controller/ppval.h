//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  ppval.h
//
//  Code generation for function 'ppval'
//


#ifndef PPVAL_H
#define PPVAL_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_controller_types.h"

// Function Declarations
extern void ppval(const double pp_breaks[400], const double pp_coefs[14364],
                  double x, double v[9]);

#endif

// End of code generation (ppval.h)
