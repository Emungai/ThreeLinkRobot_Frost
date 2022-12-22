/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 16:38:26 GMT-04:00
 */

#ifndef TAU_VCONSTRAINT_SWING_HH
#define TAU_VCONSTRAINT_SWING_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void tau_VConstraint_Swing_raw(double *p_output1, const double *var1,const double *var2);

  inline void tau_VConstraint_Swing(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);
    assert_size_matrix(var2, 2, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 1, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    tau_VConstraint_Swing_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // TAU_VCONSTRAINT_SWING_HH
