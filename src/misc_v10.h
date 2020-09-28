
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
// define very large/small numbers for catching overflow/underflow problems
const int OVERFLO_INT = INT_MAX/100;
const double OVERFLO_DOUBLE = DBL_MAX/100;
const double UNDERFLO_DOUBLE = DBL_MIN/100;

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed
// to be -inf.
double log_sum(double logA, double logB);

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
// sum boolean values and return integer
int sum_bool(const std::vector<bool> &x_vec);

//------------------------------------------------
// mean of vector
template<class TYPE>
double mean(const std::vector<TYPE> &x) {
  return sum(x)/double(x.size());
}

//------------------------------------------------
// min of vector
template<class TYPE>
TYPE min(const std::vector<TYPE> & x) {
  return *min_element(x.begin(), x.end());
}

//------------------------------------------------
// max of vector
template<class TYPE>
TYPE max(std::vector<TYPE> x) {
  return *max_element(x.begin(), x.end());
}

//------------------------------------------------
// square function
template<class TYPE>
TYPE sq(const TYPE x) {
  return x*x;
}

//------------------------------------------------
// analogue of R function seq() for integers
std::vector<int> seq_int(int from, int to, int by = 1);

//------------------------------------------------
// Euclidian distance between points in 2 dimensions
template<class TYPE>
double dist_euclid_2d(TYPE x1, TYPE y1, TYPE x2, TYPE y2) {
  return sqrt(sq(x1 - x2) + sq(y1 - y2));
}

//------------------------------------------------
// push back multiple values to vector
template<class TYPE>
void push_back_multiple(std::vector<TYPE> &lhs, const std::vector<TYPE> &rhs) {
  lhs.insert(lhs.end(), rhs.begin(), rhs.end());
}

//------------------------------------------------
// erase particular element of a vector using efficient method. Warning - does
// not preserve original order of vector
template<class TYPE>
void quick_erase(std::vector<TYPE> &v, int index) {
  v[index] = v.back();
  v.pop_back();
}

//------------------------------------------------
// erase-remove idiom, for erasing all instances of a particular value from
// container
template<class TYPE>
void erase_remove(std::vector<TYPE> &v, TYPE x) {
  v.erase(remove(v.begin(), v.end(), x), end(v));
}

//------------------------------------------------
// test whether value can be found in vector
template<class TYPE>
bool is_in_vector(TYPE x, const std::vector<TYPE> &v) {
  return find(v.begin(), v.end(), x) != v.end();
}

//------------------------------------------------
// return unique values in a vector
template<class TYPE>
std::vector<TYPE> unique(const std::vector<TYPE> &v) {
  std::vector<TYPE> ret;
  for (int i = 0; i < int(v.size()); ++i) {
    if (find(ret.begin(), ret.end(), v[i]) == ret.end()) {
      ret.push_back(v[i]);
    }
  }
  return ret;
}

//------------------------------------------------
// remove duplicated values from an already-sorted vector
template<class TYPE>
void remove_duplicates(std::vector<TYPE> &v) {
  int v_size = int(v.size());
  if (v_size == 1) {
    return;
  }
  for (int i = 1; i < v_size; ++i) {
    if (i >= int(v.size())) {
      break;
    }
    if (v[i] == v[i-1]) {
      v.erase(v.begin()+i);
      i--;
    }
  }
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

//------------------------------------------------
// read values from comma-separated text file to vector<int>
std::vector<int> file_to_vector_int(std::string file_path);

//------------------------------------------------
// read values from comma-separated text file to vector<double>
std::vector<double> file_to_vector_double(std::string file_path);

//------------------------------------------------
// read values from text file to vector<vector<double>>. Text file should be
// delimited by first line break, then comma. Lines do not all need to be same
// length, i.e. jagged matrices are allowed.
std::vector<std::vector<double>> file_to_matrix_double(std::string file_path);

//------------------------------------------------
// calculate Cholesky decomposition of positive definite matrix sigma
void cholesky(std::vector<std::vector<double>> &chol, const std::vector<std::vector<double>> &sigma);

//------------------------------------------------
// given [x,y] coordinates of points, and a series of values x_pred at which to
// return, calculate cubic spline interpolation between coordinates and save
// result into y_pred vector. Both x and x_pred must be increasing, and all
// values in x_pred must be inside (or equal to) x.
void cubic_spline(std::vector<double> &x, std::vector<double> &y,
                  std::vector<double> &x_pred, std::vector<double> &y_pred);
