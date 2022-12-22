/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 16:38:23 GMT-04:00
 */

#ifndef YA_VCONSTRAINT_SWING_HH
#define YA_VCONSTRAINT_SWING_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void ya_VConstraint_Swing_raw(double *p_output1, const double *var1);

  inline void ya_VConstraint_Swing(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 2, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    ya_VConstraint_Swing_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // YA_VCONSTRAINT_SWING_HH
