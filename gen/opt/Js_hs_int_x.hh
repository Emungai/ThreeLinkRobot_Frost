/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:31 GMT-04:00
 */

#ifndef JS_HS_INT_X_HH
#define JS_HS_INT_X_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace Pattern[ThreeLink, Blank[system]]
{

  void Js_hs_int_x_raw(double *p_output1, const double *var1);

  inline void Js_hs_int_x(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 1, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 70, 2);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    Js_hs_int_x_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // JS_HS_INT_X_HH
