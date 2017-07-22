#ifndef PARAMS
#define PARAMS

#include <iostream>
#include <fstream>
#include "Eigen/Eigen"

#include "Utils.cpp"

using namespace std;
using namespace Eigen;

struct ParamsBLY {
  double mu, sigma, r, p, h, K1, K2, c1, c2;
  struct Solver {
    double rho1, rho2;
    double alpha, beta;
    double aStar, bStar;
    double aHat, bHat;
    double aAlpha_of_bHat, bBeta_of_aHat;
    double aTemp, bTemp;
    double ATemp, BTemp;
    double aStar_lo, bStar_lo, aStar_hi, bStar_hi;
    double a_lo, b_lo, a_hi, b_hi;
  } solver;

  void print() {
    ostream& out = io.get_io();

    string ret = "";
    ret = ret
      + ", mu: " + to_string(mu)
      + ", sigma: " + to_string(sigma)
      + ", r: " + to_string(r)
      + ", p: " + to_string(p)
      + ", h: " + to_string(h)
      + ", K1: " + to_string(K1)
      + ", K2: " + to_string(K2)
      + ", c1: " + to_string(c1)
      + ", c2: " + to_string(c2)
      + "\n"
      + "Solver: "
      + ", rho1: " + to_string(solver.rho1)
      + ", rho2: " + to_string(solver.rho2)
      + ", alpha: " + to_string(solver.alpha)
      + ", beta: " + to_string(solver.beta)
      + ", aTemp: " + to_string(solver.aTemp)
      + ", ATemp: " + to_string(solver.ATemp)
      + ", BTemp: " + to_string(solver.BTemp)
      + ", bTemp: " + to_string(solver.bTemp)
      + "\n";
    out << ret << endl;
  }

  void print_solution(int iter = 0) {
    ostream& out = io.get_io();

    string ret = "";
    ret = ret
      + "(a, A, B, b)"
      + (iter == 0 ? ": " : " Iteration " + to_string(iter) + ": ")
      + to_string(solver.aTemp)
      + ", " + to_string(solver.ATemp)
      + ", " + to_string(solver.BTemp)
      + ", " + to_string(solver.bTemp)
      + "\n";
    out << ret;
  }
};

struct Params {
  double mu, sigma, beta, p, q, K, L, k, l;

  void print() {
    ostream& out = io.get_io();

    string ret = "";
    ret = ret
      + ", mu: " + to_string(mu)
      + ", sigma: " + to_string(sigma)
      + ", beta: " + to_string(beta)
      + ", p: " + to_string(p)
      + ", q: " + to_string(q)
      + ", K: " + to_string(K)
      + ", L: " + to_string(L)
      + ", k: " + to_string(k)
      + ", l: " + to_string(l);
    out << ret << endl;
  }

  ParamsBLY to_ParamsBLY() {
    ParamsBLY ret = {
        .mu = mu, .sigma = sigma, .r = beta, .p = p, .h = q, .K1 = K, .K2 = L, .c1 = k, .c2 = l
      };
    return ret;
  }
};

struct Control {
  double d, D, U, u;

  void print() {
    ostream& out = io.get_io();

    string ret = "";
    out << ret + to_string(d) + ", " + to_string(D) + ", " + to_string(U) + ", " + to_string(u)  << endl;
  }
};

struct Mesh {
  double dx;
  int N, nDofs;
};

struct Valfun {
  Mesh mesh;

  VectorXd Y;
  VectorXd dY;
};

struct Case {
  Params params;
  Control control;
};

#endif
