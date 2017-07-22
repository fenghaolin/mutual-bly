#ifndef CASES_H
#define CASES_H

#include <iostream>
#include <vector>

#include "Params.cpp"

using namespace std;

vector<Case> get_cases() {
  std::vector<Case> cases;

  // case 0
  cases.push_back({
    .params = {
      .mu = -0.2, .sigma = 0.6, .beta = 0.01, .p = 0.12, .q = 0.08, .K = 0.14, .L = 0.14, .k = 0.85, .l = 0.85
    },
    .control = {
      .d = -3.0, .D = 1.0, .U = 8.0, .u = 10
    }
  });

  // case 1
  cases.push_back({
    .params = {
      .mu = 0.02, .sigma = 0.3, .beta = 0.03, .p = 0.07, .q = 0.3, .K = 4.0, .L = 11.0, .k = 1.0, .l = 1.5
    },
    .control = {
      .d = -20.0, .D = 1.0, .U = 8.0, .u = 15
    }
  });

  // case 2
  cases.push_back({
    .params = {
      .mu = 0.04, .sigma = 0.3, .beta = 0.01, .p = 100.0, .q = 80, .K = 1.25, .L = 1.3, .k = 1.5, .l = 2.0
    },
    .control = {
      .d = -1.0, .D = -0.5, .U = 5.0, .u = 10.0
    }
  });

  // case 3
  cases.push_back({
    .params = {
      .mu = 0.02, .sigma = 0.9, .beta = 0.02, .p = 300.0, .q = 150.0, .K = 150000.0, .L = 35000.0, .k = 7000.0, .l = 6000.0
    },
    .control = {
      .d = -60.0, .D = -20.0, .U = 20.0, .u = 80.0
    }
  });

  // case 4
  cases.push_back({
    .params = {
      .mu = 0.02, .sigma = 0.9, .beta = 0.02, .p = 160.0, .q = 150.0, .K = 3600.0, .L = 35000.0, .k = 7000.0, .l = 6000
    },
    .control = {
      .d = -60.0, .D = -20.0, .U = 20.0, .u = 80.0
    }
  });

  // case 5
  cases.push_back({
    .params = {
      .mu = 0.25, .sigma = 0.55, .beta = 0.05, .p = 15, .q = 2.0, .K = 500.0, .L = 280.0, .k = 16.0, .l = 4.0
    },
    .control = {
      .d = -60.0, .D = -20.0, .U = 20.0, .u = 80.0
    }
  });

  return cases;
}

#endif
