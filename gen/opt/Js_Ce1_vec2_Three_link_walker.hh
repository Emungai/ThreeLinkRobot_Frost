/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:44 GMT-04:00
 */

#ifndef JS_CE1_VEC2_THREE_LINK_WALKER_HH
#define JS_CE1_VEC2_THREE_LINK_WALKER_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace Pattern[ThreeLink, Blank[system]]
{

  void Js_Ce1_vec2_Three_link_walker_raw(double *p_output1, const double *var1);

  inline void Js_Ce1_vec2_Three_link_walker(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 1, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 31, 2);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    Js_Ce1_vec2_Three_link_walker_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // JS_CE1_VEC2_THREE_LINK_WALKER_HH
