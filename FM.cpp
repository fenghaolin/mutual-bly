#ifndef FM_H
#define FM_H

#include <iostream>
#include "Eigen/Eigen"

#include "Solver.cpp"

const double MAX_NUM = 1.0e+100;

bool change_du(Control &control, const Valfun &V, const Params &params) {
  bool ret = false;

  int N = V.mesh.N;
  double dx = V.mesh.dx;
  double x0 = control.d;

  if ( !( V.dY(0) + params.k < 0 )) { // it could violate the condition due to last move
    for (int i = 0; i < N ; i++) {
      if (V.dY(i+1) + params.k < 0) {
        if (i != 0) {
          control.d = x0 + i * dx;
          ret = true;
        }
        break;
      }
    }
  }

  if ( !( -V.dY(N-1) + params.l < 0)) { // it could violate the condition due to last move
    for (int i = N-1; i > 0 ; i--) {
      if (-V.dY(i-1) + params.l < 0) {
        if (i != N-1) {
          control.u = x0 + i * dx;
          ret = true;
        }
        break;
      }
    }
  }

  return ret;
}

bool change_DU(Control &control, const Valfun &V, const Params &params) {
  bool ret = false;

  int N = V.mesh.N;
  double dx = V.mesh.dx;
  double x0 = control.d;

  double min_D = MAX_NUM;
  double min_U = MAX_NUM;
  double min_D_x = x0;
  double min_U_x = x0;

  double compare;
  for (int i = 0; i < N ; i++) {
    double x = x0 + dx * i;

    compare = V.Y(i) + params.k * x;
    if (min_D > compare) {
      min_D = compare;
      min_D_x = x;
    }

    compare = V.Y(i) - params.l * x;
    if (min_U > compare) {
      min_U = compare;
      min_U_x = x;
    }
  }

  if (abs(control.D - min_D_x) > dx) {
    control.D = min_D_x;
    ret = true;
  }

  if (abs(control.U - min_U_x) > dx) {
    control.U = min_U_x;
    ret = true;
  }

  if (!(min_D_x < min_U_x)) {
    throw "D and U are crossing over";
  }

  return ret;
}

void solve_control(Valfun &V, const Control &control, const Params &params, const double tolerance) {
  build_mesh(V.mesh, control, tolerance);
  solve(V, params, control);
  calc_diff_V(V);
}

void FM(Control &control, const Params &params, const double tolerance) {
  Valfun V;

  for (int i=1; true; i++) {
    bool is_goal = true;

    solve_control(V, control, params, tolerance);
    if (change_du(control, V, params)) {
      is_goal = false;
    }

    solve_control(V, control, params, tolerance);
    if (change_DU(control, V, params)) {
      is_goal = false;
    }

    if (is_goal) {
      break;
    }

    io.get_io() << "(d, D, U, u) Iteration " << i << ": ";
    control.print();
  }
}

void test_FM() {
  bool is_fail = false;

  Params params = {.mu = -0.2, .sigma = 0.6, .beta = 0.01, .p = 0.12, .q = 0.08, .K = 0.14, .L = 0.14, .k = 0.85, .l = 0.85};
  Control control = {.d = -3.0, .D = 1.0, .U = 8.0, .u = 10.0};
  double tolerance = 0.5;

  FM(control, params, tolerance);

  if (is_diff(control.d, -1.96)) is_fail = true;
  if (is_diff(control.u, 6.36)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

#endif
