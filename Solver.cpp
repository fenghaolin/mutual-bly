#ifndef SOLVER
#define SOLVER

#include <iostream>
#include "Eigen/Eigen"
#include "Eigen/SparseLU"

#include "Utils.cpp"
#include "Params.cpp"
#include "Init.cpp"

using namespace std;
using namespace Eigen;

void build_mesh(Mesh &mesh, const Control &control, const double tolerance) {
    double Length = control.u - control.d;
    mesh.N = ceil(Length / tolerance);
    mesh.nDofs = mesh.N + 1;
    mesh.dx = Length / (mesh.N-1);
}

bool solve(Valfun &V, const Params &params, const Control &control) {
    Mesh mesh = V.mesh;

    VectorXd RHS(mesh.nDofs);
    SparseMatrix<double> AMATRX(mesh.nDofs, mesh.nDofs);

    // build components of the linear system
    add_diff_rows(AMATRX, params, mesh.N, mesh.dx);
      // need to add after the diff rows
    add_first_row(AMATRX, control, mesh.dx);
    add_last_row(AMATRX, control, mesh.dx);

    init_RHS(RHS, params, control, mesh.N, mesh.dx);
    // Solve the linear system
    SparseLU<SparseMatrix<double> > solver;

    solver.compute(AMATRX);
    if(solver.info()!=Success)
    {
        cout<<"Decomposition failed!!!"<<endl;
        return false;
    }

    V.Y=solver.solve(RHS);
    if(solver.info()!=Success)
    {
        cout<<"Solving failed!!!"<<endl;
        return false;
    }

    return true;
}

void calc_diff_V(Valfun &V) {
  V.dY = VectorXd(V.mesh.N);
  for (int i = 0; i < V.mesh.N; i++) {
    V.dY(i) = ( V.Y(i+1) - V.Y(i) ) / V.mesh.dx;
  }
}

void test_solve() {
    bool is_fail = false;

    Params params = {.mu = -0.2, .sigma = 0.6, .beta = 0.01, .p = 0.12, .q = 0.08, .K = 0.14, .L = 0.14, .k = 0.85, .l = 0.85};
    Control control = {.d = -3.0, .D = 1.0, .U = 8.0, .u = 10.0};
    // double tolerance = 0.001;
    double tolerance = 0.5;
    Valfun V;
    build_mesh(V.mesh, control, tolerance);

    solve(V, params, control);
    calc_diff_V(V);

    if (is_diff(V.Y(0), 31.411657)) is_fail = true;
    if (is_diff(V.Y(26), 33.033516)) is_fail = true;

    if (is_fail) throw_error(__FUNCTION__);
}

#endif
