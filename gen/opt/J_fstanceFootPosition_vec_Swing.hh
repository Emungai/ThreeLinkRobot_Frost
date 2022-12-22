/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:58 GMT-04:00
 */

#ifndef J_FSTANCEFOOTPOSITION_VEC_SWING_HH
#define J_FSTANCEFOOTPOSITION_VEC_SWING_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace Pattern[ThreeLink, Blank[system]]
{

  void J_fstanceFootPosition_vec_Swing_raw(double *p_output1, const double *var1,const double *var2);

  inline void J_fstanceFootPosition_vec_Swing(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);
    assert_size_matrix(var2, 2, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 10, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    J_fstanceFootPosition_vec_Swing_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // J_FSTANCEFOOTPOSITION_VEC_SWING_HH
