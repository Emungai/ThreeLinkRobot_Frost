/*
 * Automatically Generated from Mathematica.
 * Fri 10 Jun 2022 15:49:23 GMT-04:00
 */

#ifndef VELOCITY_STANCE_FOOT_FUNCTION_HH
#define VELOCITY_STANCE_FOOT_FUNCTION_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymExpression
{

  void velocity_stance_foot_function_raw(double *p_output1, const double *var1,const double *var2);

  inline void velocity_stance_foot_function(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);
    assert_size_matrix(var2, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 2, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    velocity_stance_foot_function_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // VELOCITY_STANCE_FOOT_FUNCTION_HH
