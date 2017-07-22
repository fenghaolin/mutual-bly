#include <iostream>

#include "Params.cpp"
#include "Solver.cpp"
#include "FM.cpp"
#include "Benchmark.cpp"

using namespace std;

int main()
{
  io.set_mode(false);
  ostream& out = io.get_io();

  double tolerance = 0.001;
  out << "Testing FM method with tolerance: " << tolerance << endl;
  fm_performance(tolerance);

  tolerance = 1e-5;
  out << "\nTesting BLY method (SAGC initialization) with tolerance: " << tolerance << endl;
  bly_performance(tolerance);
  out << "\nTesting BLY method (without SAGC initialization) with tolerance: " << tolerance << endl;
  bly_performance(tolerance, false);

  return 0;
}

void test_all() {
  test_add_first_row();
  test_add_last_row();
  test_add_diff_rows();
  test_init_RHS();
  test_solve();
  test_FM();

  test_rho1_fun();
  test_rho2_fun();
  test_Z_fun();
  test_solve_root();
  test_phi_fun();
  test_u_fun();
  test_xi1_fun();
  test_F_fun();
  test_a_of_b_fun();
  test_G_fun();
  test_b_of_a_fun();
  test_aStar_fun();
  test_bStar_fun();
  test_solve_muti_dim_roots();
}

/*
.L libgsl.so
.L libgslcblas.so
.L Main.cpp
test_all()

std::cout << MatrixXd(AMATRX) << std::endl;
std::cout << MatrixXd(AMATRX).row(0) << std::endl;
 */
