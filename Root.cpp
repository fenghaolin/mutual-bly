#include <iostream>
#include <vector>

#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_roots.h>
#include <gsl_multiroots.h>

#include "Funs.cpp"
#include "Params.cpp"

using namespace std;

typedef double (* vFunCall)(double x, void *params);
typedef int (* vFunCall3d)(const gsl_vector * x, void *params, gsl_vector * f);

bool solve_root(double &solution, vFunCall target, ParamsBLY &p, double x_lo, double x_hi, double tolerance) {
  gsl_function F;
  F.function = target;
  F.params = &p;

  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, tolerance);

    if (status == GSL_SUCCESS) {
      solution = r;
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return status == GSL_SUCCESS;
}

void test_solve_root() {
  bool is_fail = false;

  ParamsBLY p = {
      .mu = 0.02, .sigma = 0.3, .r= 0.03, .p = 0.07, .h = 0.3, .K1 = 4.0, .K2 = 11.0, .c1 = 1.0, .c2 = 1.5
    };
  calc_rho1_rho2(p);

  double solution;

  solve_root(solution, alpha_fun, p, 0.1, 15.0, 0.001);
  if (is_diff(solution, 5.4467092)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

void print_state (size_t iter, gsl_multiroot_fsolver * s) {
  printf ("iter = %3u x = % .20f % .20f "
          "f(x) = % .3e % .3e\n",
          (int)iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}

bool solve_multi_dim_roots(vector<double> &solution, vFunCall3d target, ParamsBLY &p, vector<double> x_init, vector<double> x_lo, vector<double> x_hi, double tolerance) {
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = solution.size();
  gsl_multiroot_function f = {target, n, &p};

  gsl_vector *x = gsl_vector_alloc (n);

  for( int i=0; i<n; i++ ) {
    gsl_vector_set (x, i, x_init[i]);
  }

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      // print_state(iter, s);

      if (status == GSL_SUCCESS) {
        for(int i=0; i<solution.size(); i++) {
          solution[i] = gsl_vector_get (s->x, i);
        }
      }

      if (status)   /* check if solver is stuck */
        break;

      status =
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

  return status == GSL_SUCCESS;
}

void test_solve_muti_dim_roots() {
  bool is_fail = false;

  double tolerance = 1e-5;

  ParamsBLY p = {
      .mu = 0.02, .sigma = 0.3, .r= 0.03, .p = 0.07, .h = 0.3, .K1 = 4.0, .K2 = 11.0, .c1 = 1.0, .c2 = 1.5
    };

  calc_rho1_rho2(p);
  p.solver.alpha = 5.4467092;
  p.solver.beta = 3.499001973;

  vector<double> solution(2);
  vector<double> x_init = {-p.solver.alpha, p.solver.beta};
  vector<double> x_lo(2);
  vector<double> x_hi(2);

  bool is_succeed = solve_multi_dim_roots(solution, crossing_aAlpha_bBeta_fun, p, x_init, x_lo, x_hi, tolerance);
  if (!is_succeed) is_fail = true;
  if (is_diff(solution[0], -7.996068)) is_fail = true;
  if (is_diff(solution[1], 3.022831)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}
