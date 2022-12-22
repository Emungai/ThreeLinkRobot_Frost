/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 15:57:09 GMT-04:00
 */

#ifndef FSTANCEFOOTEND_MAP_SWING_HH
#define FSTANCEFOOTEND_MAP_SWING_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void fStanceFootEnd_map_Swing_raw(double *p_output1, const double *var1);

  inline void fStanceFootEnd_map_Swing(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 5, 3);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    fStanceFootEnd_map_Swing_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // FSTANCEFOOTEND_MAP_SWING_HH
