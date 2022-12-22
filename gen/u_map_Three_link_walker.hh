/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 16:38:28 GMT-04:00
 */

#ifndef U_MAP_THREE_LINK_WALKER_HH
#define U_MAP_THREE_LINK_WALKER_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void u_map_Three_link_walker_raw(double *p_output1, const double *var1);

  inline void u_map_Three_link_walker(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 5, 2);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    u_map_Three_link_walker_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // U_MAP_THREE_LINK_WALKER_HH
