#include <iostream>
#include "Eigen/Eigen"

#include "Utils.cpp"
#include "Params.cpp"

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;

// Initial first row of the matrix
void add_first_row(SparseMatrix<double> &AMATRX, const Control &control, double dx) {
  int n_D = round((control.D - control.d) / dx) + 1;

  AMATRX.insert(0, 1 - 1) = 1;
  AMATRX.insert(0, n_D - 1) = -1;
}

void test_add_first_row() {
  bool is_fail = false;

  Control control = {.d = 1, .D = 2, .U = 3, .u = 4};
  int N = 5;
  SparseMatrix<double> AMATRX(N+1, N+1);
  double dx = 1;

  add_first_row(AMATRX, control, dx);

  if (is_diff(AMATRX.coeff(0, 0), 1)) is_fail = true;
  if (is_diff(AMATRX.coeff(0, 1), -1)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

// Initial A_Np1 N+1'th row of the matrix
void add_last_row(SparseMatrix<double> &AMATRX, const Control &control, double dx) {
  int N = AMATRX.rows() - 1;
  int n_U = round((control.U - control.d) / dx) + 1;

  AMATRX.insert(N, n_U - 1) = -1;
  AMATRX.insert(N, N - 1) = 1;
}

void test_add_last_row() {
  bool is_fail = false;

  Control control = {.d = 1, .D = 2, .U = 3, .u = 4};
  int N = 5;
  SparseMatrix<double> AMATRX(N+1, N+1);
  double dx = 1;

  add_last_row(AMATRX, control, dx);

  if (is_diff(AMATRX.coeff(5, 2), -1)) is_fail = true;
  if (is_diff(AMATRX.coeff(5, 4), 1)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

// Initial 
void init_RHS(VectorXd &RHS, const Params &params, const Control &control, const int &N, const double &dx) {
  RHS(0) = params.K + params.k * (control.D - control.d);
  RHS(N) = params.L + params.l * (control.u - control.U);

  for (int n = 2; n <= N; n++){
    double x = control.d + (n - 1) * dx;
    RHS(n-1) =  2 * ( x < 0 ? params.p : -params.q ) * x / (params.sigma * params.sigma);
  }
}

void test_init_RHS() {
  bool is_fail = false;

  Params params = {.mu = 1, .sigma = 1, .beta = 1, .p = 1, .q = 1, .K = 1, .L = 1, .k = 1, .l = 1};
  Control control = {.d = -1, .D = 0, .U = 1, .u = 2};
  int N = 7;
  double dx = 0.5;
  VectorXd RHS(N+1);

  init_RHS(RHS, params, control, N, dx);
  // cout << RHS << endl

  if (is_diff(RHS(0), 2)) is_fail = true;
  if (is_diff(RHS(N), 2)) is_fail = true;

  if (is_diff(RHS(1), -1)) is_fail = true;
  if (is_diff(RHS(6), -4)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

// cautious: this will clear previous assigned elements in AMATRX
void add_diff_rows(SparseMatrix<double> &AMATRX, const Params &params, const int &N, const double &dx) {
  vector<T> tripletList;

  double A_n_nm1 =  1.0 / (dx * dx) - params.mu / (dx * params.sigma * params.sigma); // a_{n,n-1}
  double A_n_n =  -2.0 / (dx * dx) - 2.0 * params.beta / (params.sigma * params.sigma); // a_{n,n}
  double A_n_np1 =  1.0 / (dx * dx) + params.mu / (dx * params.sigma * params.sigma); // a_{n,n+1}

  for (int n=2; n <= N; n++) {
    tripletList.push_back(T(n - 1, n - 2 , A_n_nm1));
    tripletList.push_back(T(n - 1, n - 1 , A_n_n));
    tripletList.push_back(T(n - 1, n , A_n_np1));
  }

  AMATRX.setFromTriplets(tripletList.begin(),tripletList.end());
  tripletList.clear();
}

void test_add_diff_rows() {
  bool is_fail = false;

  Params params = {.mu = 1, .sigma = 1, .beta = 1, .p = 1, .q = 1, .K = 1, .L = 1, .k = 1, .l = 1};
  int N = 5;
  double dx = 0.1;
  VectorXd A_1(N+1), A_Np1(N+1), RHS(N+1);
  SparseMatrix<double> AMATRX(N+1, N+1);

  add_diff_rows(AMATRX, params, N, dx);

  if (is_diff(AMATRX.coeff(1, 0), 90)) is_fail = true;
  if (is_diff(AMATRX.coeff(1, 1), -202)) is_fail = true;
  if (is_diff(AMATRX.coeff(1, 2), 110)) is_fail = true;
  if (is_diff(AMATRX.coeff(1, 3), 0)) is_fail = true;

  if (is_diff(AMATRX.coeff(4, 2), 0)) is_fail = true;
  if (is_diff(AMATRX.coeff(4, 3), 90)) is_fail = true;
  if (is_diff(AMATRX.coeff(4, 4), -202)) is_fail = true;
  if (is_diff(AMATRX.coeff(4, 5), 110)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}
