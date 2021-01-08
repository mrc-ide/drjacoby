
#include "misc.h"

#include <math.h>
#include <fstream>
#include <sstream>

using namespace std;

//------------------------------------------------
// basic sum over elements in a vector
// sum
// DEFINED IN HEADER

//------------------------------------------------
// update timer and optionally print time difference
void chrono_timer(chrono::high_resolution_clock::time_point &t0, string message_before, bool print_diff) {
  
  // calculate elapsed time
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t1-t0);
  double time_double = time_span.count();
  
  // print time difference
  if (print_diff) {
    message_before += to_string(time_double) + " seconds";
    print(message_before);
  }
  
  // update timer to current time
  t0 = t1;
}

//------------------------------------------------
// helper function for printing a single value
// print
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a vector or set
// print_vector
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a matrix
// print_matrix
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a 3D array
// print_array
// DEFINED IN HEADER

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void print_stars(int n, string title) {
#ifdef RCPP_ACTIVE
  Rcpp::Rcout << title << " ";
  for (int i = 0; i < n; ++i) {
    Rcpp::Rcout << "*";
  }
  Rcpp::Rcout << "\n";
#else
  std::cout << title << " ";
  for (int i = 0; i < n; ++i) {
    std::cout << "*";
  }
  std::cout << "\n";
#endif
}

//------------------------------------------------
// print "foo", with optional number e.g. "foo2"
void foo(int n) {
#ifdef RCPP_ACTIVE
  if (n == 0) {
    Rcpp::Rcout << "foo\n";
  } else {
    Rcpp::Rcout << "foo" << n << "\n";
  }
  R_FlushConsole();
#else
  if (n == 0) {
    std::cout << "foo\n";
  } else {
    std::cout << "foo" << n << "\n";
  }
#endif
}

//------------------------------------------------
// print "bar", with optional number e.g. "bar2"
void bar(int n) {
#ifdef RCPP_ACTIVE
  if (n == 0) {
    Rcpp::Rcout << "bar\n";
  } else {
    Rcpp::Rcout << "bar" << n << "\n";
  }
  R_FlushConsole();
#else
  if (n == 0) {
    std::cout << "bar\n";
  } else {
    std::cout << "bar" << n << "\n";
  }
#endif
}

//------------------------------------------------
// print "foobar", with optional number e.g. "foobar2"
void foobar(int n) {
#ifdef RCPP_ACTIVE
  if (n == 0) {
    Rcpp::Rcout << "foobar\n";
  } else {
    Rcpp::Rcout << "foobar" << n << "\n";
  }
  R_FlushConsole();
#else
  if (n == 0) {
    std::cout << "foobar\n";
  } else {
    std::cout << "foobar" << n << "\n";
  }
#endif
}

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to bool format.
int rcpp_to_bool(SEXP x) {
  return Rcpp::as<bool>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to int format.
int rcpp_to_int(SEXP x) {
  return Rcpp::as<int>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to double format.
double rcpp_to_double(SEXP x) {
  return Rcpp::as<double>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to string format.
string rcpp_to_string(SEXP x) {
  return Rcpp::as<string>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<bool> format.
vector<bool> rcpp_to_vector_bool(SEXP x) {
  return Rcpp::as<vector<bool>>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<int> format.
vector<int> rcpp_to_vector_int(SEXP x) {
  return Rcpp::as<vector<int>>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<double> format.
vector<double> rcpp_to_vector_double(SEXP x) {
  return Rcpp::as<vector<double>>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<string> format.
vector<string> rcpp_to_vector_string(SEXP x) {
  return Rcpp::as<vector<string>>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<bool>> format.
vector<vector<bool>> rcpp_to_matrix_bool(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<bool> > x_mat(nrow);
  for (int i = 0; i < nrow; ++i) {
    x_mat[i] = Rcpp::as<vector<bool>>(x[i]);
  }
  return x_mat;
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<int>> format.
vector<vector<int>> rcpp_to_matrix_int(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<int> > x_mat(nrow);
  for (int i = 0; i < nrow; ++i) {
    x_mat[i] = Rcpp::as<vector<int>>(x[i]);
  }
  return x_mat;
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<double>> format.
vector<vector<double>> rcpp_to_matrix_double(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<double> > x_mat(nrow);
  for (int i = 0; i < nrow; ++i) {
    x_mat[i] = Rcpp::as<vector<double>>(x[i]);
  }
  return x_mat;
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<string>> format.
vector<vector<string>> rcpp_to_matrix_string(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<string> > x_mat(nrow);
  for (int i = 0; i < nrow; ++i) {
    x_mat[i] = Rcpp::as<vector<string>>(x[i]);
  }
  return x_mat;
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
vector<vector<vector<int>>> rcpp_to_array_int(Rcpp::List x) {
  int n1 = int(x.size());
  vector<vector<vector<int>>> ret(n1);
  for (int i = 0; i < n1; ++i) {
    Rcpp::List x_i = x[i];
    int n2 = int(x_i.size());
    ret[i] = vector<vector<int>>(n2);
    for (int j = 0; j < n2; ++j) {
      ret[i][j] = Rcpp::as<vector<int>>(x_i[j]);
    }
  }
  return ret;
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
vector<vector<vector<double>>> rcpp_to_array_double(Rcpp::List x) {
  int n1 = int(x.size());
  vector<vector<vector<double>>> ret(n1);
  for (int i = 0; i < n1; ++i) {
    Rcpp::List x_i = x[i];
    int n2 = int(x_i.size());
    ret[i] = vector<vector<double>>(n2);
    for (int j = 0; j < n2; ++j) {
      ret[i][j] = Rcpp::as<vector<double>>(x_i[j]);
    }
  }
  return ret;
}
#endif
