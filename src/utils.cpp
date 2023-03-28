#include <cpp11.hpp>
#include <vector>

double sum(std::vector<double> x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}
