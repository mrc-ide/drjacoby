#include <cpp11.hpp>
#include <vector>
#include <chrono>

double sum(const std::vector<double> x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}

bool value_check(const double x){
  bool is_nan = std::isnan(x);
  bool is_infinite = !std::isfinite(x);
  bool is_positive = x > 0;
  bool fail = is_nan || (is_infinite && is_positive);
  return fail;
}

double chrono_timer(std::chrono::high_resolution_clock::time_point &t0, std::string message_before, bool print_diff) {
  
  // calculate elapsed time
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1-t0);
  double time_double = time_span.count();
  
  // print time difference
  if (print_diff) {
    message_before += std::to_string(time_double) + " seconds\n";
    cpp11::message(message_before);
  }
  
  // update timer to current time
  t0 = t1;
  
  // return time diff
  return time_double;
}
