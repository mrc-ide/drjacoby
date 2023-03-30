#ifndef UTILS_H
#define UTILS_H

#include <cpp11.hpp>
#include <vector>
#include <chrono>

double sum(std::vector<double> x);
double chrono_timer(std::chrono::high_resolution_clock::time_point &t0, std::string message_before, bool print_diff); 
  
#endif
