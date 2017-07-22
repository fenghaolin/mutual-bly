#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <iostream>
#include <vector>

#include "Params.cpp"
#include "Solver.cpp"
#include "Cases.cpp"
#include "FM.cpp"
#include "BLY.cpp"

using namespace std;

void fm_performance(double tolerance) {
  vector<Case> cases = get_cases();
  ostream& out = io.get_io();

  for (int i = 0; i < cases.size(); i++) {
    out << endl << "Case " << i << ":" << endl;
    cases[i].params.print();

    timer.start();
    FM(cases[i].control, cases[i].params, tolerance);
    double time_spent = timer.stop();
    out << "Time spent: " << time_spent << " seconds." << endl;
  }
}

void set_init_a_by_case(ParamsBLY &p, const int i) {
  switch(i) {
    case 0: p.solver.aTemp = -3; break;
    case 1: p.solver.aTemp = -20; break;
    case 2: p.solver.aTemp = -1; break;
    case 3: p.solver.aTemp = -60; break;
    case 4: p.solver.aTemp = -60; break;
    case 5: p.solver.aTemp = -60; break;
  }
}

void bly_performance(double tolerance, bool use_SAGC = true) {
  vector<Case> cases = get_cases();
  ostream& out = io.get_io();
  ParamsBLY p;

  for (int i = 0; i < cases.size(); i++) {
    out << endl << "Case " << i << ":" << endl;
    cases[i].params.print();
    p = cases[i].params.to_ParamsBLY();

    if (!use_SAGC) {
      set_init_a_by_case(p, i);
    }

    timer.start();
    BLY(p, tolerance, use_SAGC);
    double time_spent = timer.stop();
    out << "Time spent: " << time_spent << " seconds." << endl;
  }
}

#endif
