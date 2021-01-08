
#pragma once

// comment out this definition to switch from Rcpp to C++
#define RCPP_ACTIVE

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <chrono>
#include <cfloat>

//------------------------------------------------
// basic sum over elements in a vector
template<class TYPE>
TYPE sum(const std::vector<TYPE> &x_vec) {
  TYPE ret = 0;
  for (const auto & x : x_vec) {
    ret += x;
  }
  return ret;
}

//------------------------------------------------
// update timer and optionally print time difference
void chrono_timer(std::chrono::high_resolution_clock::time_point &t0,
                  std::string message_before = "completed in ",
                  bool print_diff = true);

//------------------------------------------------
// helper function for printing a single value or series of values
template<typename TYPE>
void print(TYPE x) {
#ifdef RCPP_ACTIVE
  Rcpp::Rcout << x << "\n";
#else
  std::cout << x << "\n";
#endif
}

template<typename TYPE, typename... Args>
void print(TYPE first, Args... args) {
#ifdef RCPP_ACTIVE
  Rcpp::Rcout << first << " ";
#else
  std::cout << first << " ";
#endif
  print(args...);
}

//------------------------------------------------
// helper function for printing contents of a vector or set
template<class TYPE>
void print_vector(const TYPE &x) {
#ifdef RCPP_ACTIVE
  for (auto i : x) {
    Rcpp::Rcout << i << " ";
  }
  Rcpp::Rcout << "\n";
#else
  for (auto i : x) {
    std::cout << i << " ";
  }
  std::cout << "\n";
#endif
}

//------------------------------------------------
// helper function for printing contents of a matrix
template<class TYPE>
void print_matrix(const std::vector<std::vector<TYPE>> &x) {
  for (int i = 0; i < x.size(); ++i) {
    print_vector(x[i]);
  }
#ifdef RCPP_ACTIVE
  Rcpp::Rcout << "\n";
#else
  std::cout << "\n";
#endif
}

//------------------------------------------------
// helper function for printing contents of a 3D array
template<class TYPE>
void print_array(const std::vector<std::vector<std::vector<TYPE>>> &x) {
#ifdef RCPP_ACTIVE
  for (int i = 0; i < x.size(); ++i) {
    Rcpp::Rcout << "--- slice " << i+1 << " ---\n";
    print_matrix(x[i]);
  }
  Rcpp::Rcout << "\n";
#else
  for (int i = 0; i < x.size(); ++i) {
    std::cout << "--- slice " << i+1 << " ---\n";
    print_matrix(x[i]);
  }
  std::cout << "\n";
#endif
}

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void print_stars(int n = 50, std::string title = "");

//------------------------------------------------
// print "foo", with optional number e.g. "foo2"
void foo(int n = 0);

//------------------------------------------------
// print "bar", with optional number e.g. "bar2"
void bar(int n = 0);

//------------------------------------------------
// print "foobar", with optional number e.g. "foobar2"
void foobar(int n = 0);

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to bool format.
int rcpp_to_bool(SEXP x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to int format.
int rcpp_to_int(SEXP x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to double format.
double rcpp_to_double(SEXP x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to string format.
std::string rcpp_to_string(SEXP x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<bool> format.
std::vector<bool> rcpp_to_vector_bool(SEXP x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<int> format.
std::vector<int> rcpp_to_vector_int(SEXP x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<double> format.
std::vector<double> rcpp_to_vector_double(SEXP x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<string> format.
std::vector<std::string> rcpp_to_vector_string(SEXP x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<bool>> format.
std::vector<std::vector<bool>> rcpp_to_matrix_bool(Rcpp::List x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<int>> format.
std::vector<std::vector<int>> rcpp_to_matrix_int(Rcpp::List x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<double>> format.
std::vector<std::vector<double>> rcpp_to_matrix_double(Rcpp::List x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<string>> format.
std::vector<std::vector<std::string>> rcpp_to_matrix_sting(Rcpp::List x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
std::vector<std::vector<std::vector<int>>> rcpp_to_array_int(Rcpp::List x);
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
std::vector<std::vector<std::vector<double>>> rcpp_to_array_double(Rcpp::List x);
#endif
