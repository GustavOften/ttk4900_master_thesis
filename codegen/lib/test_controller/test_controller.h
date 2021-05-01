//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  test_controller.h
//
//  Code generation for function 'test_controller'
//


#ifndef TEST_CONTROLLER_H
#define TEST_CONTROLLER_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "test_controller_types.h"

// Function Declarations
extern void spline_phi_not_empty_init();
extern void test_controller(double theta, double varphi, double dtheta, double
  dvarphi, double t, double luenberger_gain, double *u, double epsilon[3]);

#endif

// End of code generation (test_controller.h)
