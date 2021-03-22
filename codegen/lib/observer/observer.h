//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  observer.h
//
//  Code generation for function 'observer'
//


#ifndef OBSERVER_H
#define OBSERVER_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "observer_types.h"

// Function Declarations
extern void observer(double t, double theta, double varphi, double u, double q[2],
                     double dq[2]);
extern void observer_init();
extern void spline_phi_not_empty_init();

#endif

// End of code generation (observer.h)
