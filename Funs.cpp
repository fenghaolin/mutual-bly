#ifndef FUNS
#define FUNS

#include <gsl_math.h>
#include <gsl_multiroots.h>
#include "Utils.cpp"
#include "Params.cpp"
#include "Cases.cpp"

double rho1_fun( const ParamsBLY &p ) {
  return 
    - (p.mu + sqrt( p.mu*p.mu + 2*p.r*p.sigma*p.sigma ) ) / ( p.sigma*p.sigma );
}

void test_rho1_fun() {
  bool is_fail = false;
  ParamsBLY p = {.mu = 0.02, .r = 0.03, .sigma =0.3};
  if (is_diff(rho1_fun(p), -1.068419)) is_fail = true;

  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  if (is_diff(rho1_fun(p), -4.0)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double rho2_fun( const ParamsBLY &p ) {
  return
    - (p.mu - sqrt( p.mu*p.mu + 2*p.r*p.sigma*p.sigma ) ) / ( p.sigma*p.sigma );
}

void test_rho2_fun() {
  bool is_fail = false;
  double mu, r, sigma;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  if (is_diff(rho2_fun(p), 2.0)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

void calc_rho1_rho2(ParamsBLY &p) {
  p.solver.rho1 = rho1_fun(p);
  p.solver.rho2 = rho2_fun(p);
}

double Z_fun( const double x, const ParamsBLY &p ) {
  return
    -x + (p.solver.rho1-p.solver.rho2) / (p.solver.rho1*p.solver.rho2)
    * (exp( p.solver.rho1*x ) - 1) * (exp( p.solver.rho2*x ) - 1)
    / ( exp( p.solver.rho1*x ) - exp( p.solver.rho2*x ) );
}

void test_Z_fun() {
  bool is_fail = false;
  double x;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  x = 1.0;
  calc_rho1_rho2(p);
  if (is_diff(Z_fun(x, p), -0.361797)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double phi_fun( const double x, const ParamsBLY &p ) {
  return
    ( p.solver.rho2 * (exp( p.solver.rho1*x ) - 1) - p.solver.rho1 * (exp( p.solver.rho2*x ) - 1) )
    / ( exp( p.solver.rho2*x ) - exp( p.solver.rho1*x ) );
}

void test_phi_fun() {
  bool is_fail = false;
  double x;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  x = 1.0;
  calc_rho1_rho2(p);
  if (is_diff(phi_fun(x, p), 3.200880)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double u_fun( const double x, const ParamsBLY &p ) {
  return
    ( p.c1 + p.c2 ) * p.r * ( p.solver.rho2 - p.solver.rho1 )
    + (p.p +p.h) * ( p.solver.rho2 * (exp( p.solver.rho1*x ) - 1) - p.solver.rho1 * (exp( p.solver.rho2*x ) - 1) );
}

void test_u_fun() {
  bool is_fail = false;
  double x;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  p.c1 = 0.5, p.c2 = 0.5;
  p.p = 0.5, p.h = 0.5;
  x = 1.0;
  calc_rho1_rho2(p);
  if (is_diff(u_fun(x, p), 47.592856)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double z_fun( const double x, const ParamsBLY &p ) {
  return
    ( p.p - p.r * p.c1 )
    * (
        ( phi_fun(p.solver.alpha, p) + p.solver.rho2 ) * (exp( p.solver.rho1*x ) - 1)
        - ( phi_fun(p.solver.alpha, p) + p.solver.rho1 ) * (exp( p.solver.rho2*x ) - 1)
    );
}

void test_z_fun() {
  bool is_fail = false;
  double x;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  p.c1 = 0.5, p.c2 = 0.5;
  p.p = 0.5, p.h = 0.5;
  x = 1.0;
  p.solver.alpha = 1.0;
  calc_rho1_rho2(p);

  if (is_diff(z_fun(x, p), 47.592856)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double xi1_fun( const double x, const ParamsBLY &p ) {
  return
    ( p.h - p.r * p.c2 )
    * (
        ( phi_fun(-p.solver.beta, p) + p.solver.rho2 ) * (exp( p.solver.rho1*x ) - 1)
        - ( phi_fun(-p.solver.beta, p) + p.solver.rho1 ) * (exp( p.solver.rho2*x ) - 1)
    );
}

void test_xi1_fun() {
  bool is_fail = false;
  double x;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  p.c1 = 0.5, p.c2 = 0.5;
  p.p = 0.5, p.h = 0.5;
  x = 1.0;
  p.solver.beta = 1.0;
  calc_rho1_rho2(p);

  if (is_diff(xi1_fun(x, p), -56.448328)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double F_fun( const double a, const double A, const ParamsBLY &p ) {
  return
    p.K1 + ( ( p.p - p.r * p.c1 ) * a + ( p.h + p.r * p.c1 ) * A ) / p.r
    + (
        (p.solver.rho1 - p.solver.rho2) * ( p.p - p.r * p.c1 ) * (exp( p.solver.rho1*(A-a) ) - 1) * (exp( p.solver.rho2*(A-a) ) - 1)
        + ( p.p + p.h )
        * ( p.solver.rho1 * (exp( p.solver.rho1*(A-a) ) - 1) * ( 1 - exp( p.solver.rho2*A) )
            - p.solver.rho2 * (exp( p.solver.rho2*(A-a) ) - 1) * ( 1 - exp( p.solver.rho1*A) ) )
    )
    / ( p.r * p.solver.rho1 * p.solver.rho2 * (exp( p.solver.rho1*(A-a) ) - exp( p.solver.rho2*(A-a) ) ) );
}

void test_F_fun() {
  bool is_fail = false;
  double a, A;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  p.K1 = 10.0, p.K2 = 10.0;
  p.c1 = 0.5, p.c2 = 0.5;
  p.p = 0.5, p.h = 0.5;
  a = 1.0, A = 2.0;
  calc_rho1_rho2(p);

  if (is_diff(F_fun(a, A, p), 9.689196)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double a_of_b_fun( const double a, const double b, const double A, const ParamsBLY &p ) {
  return
    ( p.c1 + p.c2 ) * (exp( p.solver.rho1*(A-a) ) - exp( p.solver.rho2*(A-a) ) )
    - ( p.h + p.r * p.c1 ) / p.r
    * (exp( p.solver.rho1*(A-a) ) - exp( p.solver.rho2*(A-a) ) - exp( p.solver.rho1*(b-a) ) + exp( p.solver.rho2*(b-a) ) )
    + ( p.p + p.h ) / ( p.r * ( p.solver.rho2 - p.solver.rho1 ) )
    * (exp( p.solver.rho2*A + p.solver.rho1*b ) - exp( p.solver.rho1*A + p.solver.rho2*b ) )
    * ( p.solver.rho1 * exp( -p.solver.rho1*a ) - p.solver.rho2 * exp( -p.solver.rho2*a ) )
    + ( p.p - p.r * p.c1 ) / p.r
    * exp( -(p.solver.rho1+p.solver.rho2)*a )
    * (exp( p.solver.rho2*A + p.solver.rho1*b ) - exp( p.solver.rho1*A + p.solver.rho2*b ) );
}

void test_a_of_b_fun() {
  bool is_fail = false;
  double a, A, B, b;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  p.K1 = 10.0, p.K2 = 10.0;
  p.c1 = 0.5, p.c2 = 0.5;
  p.p = 0.5, p.h = 0.5;
  a = 1.0, A = 2.0;
  b = 4.0, B = 3.0;
  calc_rho1_rho2(p);

  if (is_diff(a_of_b_fun(a, b, A, p), -243.025227)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double G_fun( const double b, const double B, const ParamsBLY &p ) {
  return
    -p.K2 + ( ( p.h - p.r * p.c2 ) * b + ( p.p + p.r * p.c2 ) * B ) / p.r
    + (
        (p.solver.rho1 - p.solver.rho2) * ( p.h - p.r * p.c2 ) * (exp( p.solver.rho1*(B-b) ) - 1) * (exp( p.solver.rho2*(B-b) ) - 1)
        + ( p.p + p.h )
        * ( p.solver.rho1 * (exp( p.solver.rho1*(B-b) ) - 1) * ( 1 - exp( p.solver.rho2*B) )
            - p.solver.rho2 * (exp( p.solver.rho2*(B-b) ) - 1) * ( 1 - exp( p.solver.rho1*B) ) )
    )
    / ( p.r * p.solver.rho1 * p.solver.rho2 * (exp( p.solver.rho1*(B-b) ) - exp( p.solver.rho2*(B-b) ) ) );
}

void test_G_fun() {
  bool is_fail = false;
  double b, B;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  p.K1 = 10.0, p.K2 = 10.0;
  p.c1 = 0.5, p.c2 = 0.5;
  p.p = 0.5, p.h = 0.5;
  b = 4.0, B = 3.0;

  G_fun(b, B, p);

  if (is_diff(G_fun(b, B, p), -58.891633)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double b_of_a_fun( const double a, const double b, const double B, const ParamsBLY &p ) {
  return
    -( p.c1 + p.c2 ) * (exp( p.solver.rho1*(B-b) ) - exp( p.solver.rho2*(B-b) ) )
    + ( p.p + p.r * p.c2 ) / p.r
    * (exp( p.solver.rho1*(B-b) ) - exp( p.solver.rho2*(B-b) ) - exp( p.solver.rho1*(a-b) ) + exp( p.solver.rho2*(a-b) ) )
    + ( p.p + p.h ) / ( p.r * ( p.solver.rho1 - p.solver.rho2 ) )
    * (exp( p.solver.rho1*B + p.solver.rho2*a ) - exp( p.solver.rho2*B + p.solver.rho1*a ) )
    * ( p.solver.rho2 * exp( -p.solver.rho2*b ) - p.solver.rho1 * exp( -p.solver.rho1*b ) )
    + ( p.h - p.r * p.c2 ) / p.r
    * exp( -(p.solver.rho1+p.solver.rho2)*b )
    * (exp( p.solver.rho1*B + p.solver.rho2*a ) - exp( p.solver.rho2*B + p.solver.rho1*a ) );
}

void test_b_of_a_fun() {
  bool is_fail = false;
  double a, A, B, b;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  p.K1 = 10.0, p.K2 = 10.0;
  p.c1 = 0.5, p.c2 = 0.5;
  p.p = 0.5, p.h = 0.5;
  a = 1.0, A = 2.0;
  b = 4.0, B = 3.0;
  calc_rho1_rho2(p);

  if (is_diff(b_of_a_fun(a, b, B, p), 10849778.658909)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double aAlpha_of_b_fun( const double a, const double b, const ParamsBLY &p ) {
  return
    u_fun(b, p) - z_fun((b-a), p);
}

double bBeta_of_a_fun( const double a, const double b, const ParamsBLY &p ) {
  return
    u_fun(a, p) - xi1_fun((a-b), p);
}

double aStar_fun( const double a, const ParamsBLY &p ) {
  return
    ( p.p + p.h - (p.p - p.r*p.c1) * exp( -p.solver.rho1 * a ) )
    * pow( p.p + p.h - (p.p - p.r*p.c1) * exp( -p.solver.rho2 * a ), -p.solver.rho1 / p.solver.rho2 )
    - pow( p.h - p.r*p.c2, 1 - p.solver.rho1 / p.solver.rho2 )
    * ( p.solver.rho2 + phi_fun( -p.solver.beta, p ) ) / p.solver.rho2
    * pow( ( p.solver.rho1 + phi_fun( -p.solver.beta, p ) ) / p.solver.rho1, -p.solver.rho1 / p.solver.rho2 );
}

void test_aStar_fun() {
  bool is_fail = false;
  double x;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  p.c1 = 0.5, p.c2 = 0.5;
  p.p = 0.5, p.h = 0.5;
  x = 1.0;
  p.solver.beta = 1.0;
  calc_rho1_rho2(p);

  if (is_diff(aStar_fun(x, p), 120.320481)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

double bStar_fun( const double b, const ParamsBLY &p ) {
  return
    ( p.p + p.h - (p.h - p.r*p.c2) * exp( -p.solver.rho1 * b ) )
    * pow( p.p + p.h - (p.h - p.r*p.c2) * exp( -p.solver.rho2 * b ), -p.solver.rho1 / p.solver.rho2 )
    - pow( p.p - p.r*p.c1, 1 - p.solver.rho1 / p.solver.rho2 )
    * ( p.solver.rho2 + phi_fun( p.solver.alpha, p ) ) / p.solver.rho2
    * pow( ( p.solver.rho1 + phi_fun( p.solver.alpha, p ) ) / p.solver.rho1, -p.solver.rho1 / p.solver.rho2 );
}

void test_bStar_fun() {
  bool is_fail = false;
  double x;

  ParamsBLY p;
  p.mu = 1.0, p.r = 4.0, p.sigma =1.0;
  p.c1 = 0.5, p.c2 = 0.5;
  p.p = 0.5, p.h = 0.5;
  x = 1.0;
  p.solver.alpha = 1.0;
  calc_rho1_rho2(p);

  if (is_diff(bStar_fun(x, p), 120.320481)) is_fail = true;

  if (is_fail) throw_error(__FUNCTION__);
}

// functions for solver
double alpha_fun(double x, void *params) {
  ParamsBLY p = *( static_cast<ParamsBLY *>(params) );

  return
    Z_fun(x, p) + p.K1 * p.r / ( p.p - p.r * p.c1 );
}

double beta_fun(double x, void *params) {
  ParamsBLY p = *( static_cast<ParamsBLY *>(params) );

  return
    Z_fun(x, p) + p.K2 * p.r / ( p.h - p.r * p.c2 );
}

double solve_aStar_fun(double a, void *params) {
  ParamsBLY p = *( static_cast<ParamsBLY *>(params) );

  return
    aStar_fun(a, p);
}

double solve_bStar_fun(double b, void *params) {
  ParamsBLY p = *( static_cast<ParamsBLY *>(params) );

  return
    bStar_fun(b, p);
}

double solve_aAlpha_of_b_fun(double a, void *params) {
  ParamsBLY p = *( static_cast<ParamsBLY *>(params) );

  return
    aAlpha_of_b_fun(a, p.solver.bTemp, p);
}

double solve_bBeta_of_a_fun(double b, void *params) {
  ParamsBLY p = *( static_cast<ParamsBLY *>(params) );

  return
    bBeta_of_a_fun(p.solver.aTemp, b, p);
}

int solve_a_of_b_fun( const gsl_vector * x, void *params, gsl_vector * f ) {
  ParamsBLY p = *( static_cast<ParamsBLY *>(params) );

  const double a = gsl_vector_get (x, 0);
  const double A = gsl_vector_get (x, 1);

  const double y0 = F_fun(a, A, p);
  const double y1 = a_of_b_fun(a, p.solver.bTemp, A, p);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}

int solve_b_of_a_fun( const gsl_vector * x, void *params, gsl_vector * f ) {
  ParamsBLY p = *( static_cast<ParamsBLY *>(params) );

  const double b = gsl_vector_get (x, 0);
  const double B = gsl_vector_get (x, 1);

  const double y0 = G_fun(b, B, p);
  const double y1 = b_of_a_fun(p.solver.aTemp, b, B, p);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}

int crossing_aAlpha_bBeta_fun( const gsl_vector * x, void *params, gsl_vector * f ) {
  ParamsBLY p = *( static_cast<ParamsBLY *>(params) );

  const double a = gsl_vector_get (x, 0);
  const double b = gsl_vector_get (x, 1);

  const double y0 = aAlpha_of_b_fun(a, b, p);
  const double y1 = bBeta_of_a_fun(a, b, p);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}

#endif
