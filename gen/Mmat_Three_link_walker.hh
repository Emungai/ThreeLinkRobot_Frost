/*
 * Automatically Generated from Mathematica.
 * Wed 15 Jul 2020 15:09:02 GMT-04:00
 */

#ifndef MMAT_THREE_LINK_WALKER_HH
#define MMAT_THREE_LINK_WALKER_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void Mmat_Three_link_walker_raw(double *p_output1, const double *var1);

  inline void Mmat_Three_link_walker(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 5, 5);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    Mmat_Three_link_walker_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // MMAT_THREE_LINK_WALKER_HH
