/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:26 GMT-04:00
 */

#ifndef J_PSTANCEFOOTPOSITIONCONT_SWING_HH
#define J_PSTANCEFOOTPOSITIONCONT_SWING_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace Pattern[ThreeLink, Blank[system]]
{

  void J_pstanceFootPositionCont_Swing_raw(double *p_output1, const double *var1,const double *var2);

  inline void J_pstanceFootPositionCont_Swing(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 1, 2);
    assert_size_matrix(var2, 1, 2);

	
    // - Outputs
    assert_size_matrix(p_output1, 4, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    J_pstanceFootPositionCont_Swing_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // J_PSTANCEFOOTPOSITIONCONT_SWING_HH
