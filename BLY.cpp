#ifndef BLY_H
#define BLY_H

#define BOUNDARY 100

#include <iostream>

#include "Params.cpp"
#include "Root.cpp"
#include "Cases.cpp"

using namespace std;

void calc_alpha_beta(ParamsBLY &p, const double &tolerance) {
  double alpha, beta;
  solve_root(alpha, alpha_fun, p, 0.1, 50.0, tolerance);
  solve_root(beta, beta_fun, p, 0.1, 50.0, tolerance);
  p.solver.alpha = alpha;
  p.solver.beta = beta;
}

void calc_aStar_bStar(ParamsBLY &p, const double &tolerance) {
  solve_root(p.solver.aStar, solve_aStar_fun, p, p.solver.aStar_lo, p.solver.aStar_hi, tolerance);
  solve_root(p.solver.bStar, solve_bStar_fun, p, p.solver.bStar_lo, p.solver.bStar_hi, tolerance);
}

void calc_aHat_bHat(ParamsBLY &p, const double &tolerance) {
  vector<double> aHat_bHat(2);
  vector<double> aHat_bHat_init = {-p.solver.alpha * 1.5, p.solver.beta * 1.5};
  vector<double> aHat_bHat_lo(2);
  vector<double> aHat_bHat_hi(2);

  bool is_succeed = solve_multi_dim_roots(aHat_bHat, crossing_aAlpha_bBeta_fun, p, aHat_bHat_init, aHat_bHat_lo, aHat_bHat_hi, tolerance);
  p.solver.aHat = aHat_bHat[0];
  p.solver.bHat = aHat_bHat[1];
}

void calc_hats_values(ParamsBLY &p, const double &tolerance) {
  p.solver.aTemp = p.solver.aHat;
  p.solver.bTemp = p.solver.bHat;

  solve_root(p.solver.aAlpha_of_bHat, solve_aAlpha_of_b_fun, p, p.solver.a_lo, p.solver.a_hi, tolerance);
  solve_root(p.solver.bBeta_of_aHat, solve_bBeta_of_a_fun, p, p.solver.b_lo, p.solver.b_hi, tolerance);
}

double calc_aAlpha_of_b(ParamsBLY &p, const double &tolerance) {
  double ret;
  solve_root(ret, solve_aAlpha_of_b_fun, p, p.solver.a_lo, p.solver.a_hi, tolerance);
  return ret;
}

double calc_bBeta_of_a(ParamsBLY &p, const double &tolerance) {
  double ret;
  solve_root(ret, solve_bBeta_of_a_fun, p, p.solver.b_lo, p.solver.b_hi, tolerance);
  return ret;
}

vector<double> calc_a_of_b(const double b, ParamsBLY &p, const double &tolerance) {
  p.solver.bTemp = b;

  vector<double> a_A(2);
  vector<double> a_A_init = {p.solver.aTemp, p.solver.ATemp};
  vector<double> a_A_lo(2);
  vector<double> a_A_hi(2);

  bool is_succeed = solve_multi_dim_roots(a_A, solve_a_of_b_fun, p, a_A_init, a_A_lo, a_A_hi, tolerance);
  if (!is_succeed) {
    cout << "Fail to find the root" << endl;
  }
  return a_A;
}

vector<double> calc_b_of_a(const double a, ParamsBLY &p, const double &tolerance) {
  p.solver.aTemp = a;

  vector<double> b_B(2);
  vector<double> b_B_init = {p.solver.bTemp, p.solver.BTemp};
  vector<double> b_B_lo(2);
  vector<double> b_B_hi(2);

  bool is_succeed = solve_multi_dim_roots(b_B, solve_b_of_a_fun, p, b_B_init, b_B_lo, b_B_hi, tolerance);
  if (!is_succeed) {
    cout << "Fail to find the root" << endl;
  }
  return b_B;
}

void SAGC(ParamsBLY &p, const double &tolerance) {
  calc_aHat_bHat(p, tolerance);

  p.solver.aTemp = p.solver.aHat;
  p.solver.bTemp = p.solver.bHat;
  p.solver.ATemp = p.solver.aTemp + p.solver.alpha;
  p.solver.BTemp = p.solver.bTemp - p.solver.beta;
}

void BLY(ParamsBLY &p, const double &tolerance, bool use_SAGC = true) {
  int iter = 0;

  p.solver.a_lo = -BOUNDARY;
  p.solver.a_hi = -0.001;
  p.solver.b_lo = 0.001;
  p.solver.b_hi = BOUNDARY;

  calc_rho1_rho2(p);
  calc_alpha_beta(p, tolerance);

  if (use_SAGC) {
    SAGC(p, tolerance);
  }

  if (p.solver.ATemp < 0 && p.solver.BTemp > 0) {
    iter++;
    p.print_solution(iter);
    return;
  }

  bool is_a_converged, is_b_converged;
  double a_old, b_old;
  double aAlpha, bBeta;
  vector<double> a_A, b_B;
  do {
    is_a_converged = false;
    is_b_converged = false;
    a_old = p.solver.aTemp;
    b_old = p.solver.bTemp;

    aAlpha = calc_aAlpha_of_b(p, tolerance);
    p.solver.aTemp = aAlpha;
    p.solver.ATemp = aAlpha + p.solver.alpha;
    if (p.solver.ATemp > 0) {
      vector<double> a_A = calc_a_of_b(p.solver.bTemp, p, tolerance);
      p.solver.aTemp = a_A[0];
      p.solver.ATemp = a_A[1];
    }
    if (abs(p.solver.aTemp - a_old) < tolerance) {
      is_a_converged = true;
    }

    bBeta = calc_bBeta_of_a(p, tolerance);
    p.solver.bTemp = bBeta;
    p.solver.BTemp = bBeta - p.solver.beta;
    if (p.solver.BTemp < 0) {
      vector<double> b_B = calc_b_of_a(p.solver.aTemp, p, tolerance);
      p.solver.bTemp = b_B[0];
      p.solver.BTemp = b_B[1];
    }
    if (abs(p.solver.bTemp - b_old) < tolerance) {
      is_b_converged = true;
    }

    iter++;
    p.print_solution(iter);
  } while (!is_a_converged || !is_b_converged);
}

#endif
/*
.L libgsl.so
.L libgslcblas.so
.L BLY.cpp
 */
