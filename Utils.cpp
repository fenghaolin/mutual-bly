#ifndef UTILS
#define UTILS

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void throw_error(string message) {
  cout << "TEST Failure: " << message << endl;
  throw 1;
}

bool is_diff(double v1, double v2) {
  if ( abs(v1 - v2) > 0.000001 ) {
    cout << to_string(v1) + ", " + to_string(v2) << endl << flush;
    return true;
  }
  return false;
}

struct IO {
  string file_path = "output/benchmark.txt";
  ostream& console_out = cout;
  ofstream file_out;
  bool is_development = true;

  IO() {
    file_out.open(file_path, ios_base::app);
  }

  ostream& get_io() {
    if (is_development) {
      return console_out;
    } else {
      return file_out;
    }
  }

  void set_mode(bool is_dev){
    is_development = is_dev;
    if (!is_dev) {
      // clear the file
      file_out.close();
      file_out.open( file_path, ios::out | ios::trunc  );
      file_out.close();

      // reload the file
      file_out.open(file_path, ios_base::app);
    }
  }
};
static IO io;

struct Timer {
  double start_time;

  void start() {
    start_time = clock();
  }

  double stop() {
    return ( clock() - start_time  ) / (double) CLOCKS_PER_SEC;
  }
};
static Timer timer;

#endif
