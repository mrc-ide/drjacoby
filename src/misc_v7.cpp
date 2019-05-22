
#include "misc_v7.h"

#include <math.h>
#include <fstream>
#include <sstream>

using namespace std;

//------------------------------------------------
// define very large/small numbers for catching overflow/underflow problems
// OVERFLO_INT
// OVERFLO_DOUBLE
// UNDERFLO_DOUBLE
// DEFINED IN HEADER

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed
// to be -inf.
double log_sum(double logA, double logB) {
  if (logA-logB > 100) {
    return(logA);
  } else if (logB-logA > 100) {
    return(logB);
  }
  return (logA<logB) ? logB + log(1 + exp(logA-logB)) : logA + log(1 + exp(logB-logA));
}

//------------------------------------------------
// basic sum over elements in a vector
// sum
// DEFINED IN HEADER

//------------------------------------------------
// sum boolean values and return integer
int sum_bool(const vector<bool> &x_vec) {
  int ret = 0;
  for (int i=0; i<int(x_vec.size()); ++i) {
    ret += x_vec[i];
  }
  return ret;
}

//------------------------------------------------
// mean of vector
// mean
// DEFINED IN HEADER

//------------------------------------------------
// min of vector
// min
// DEFINED IN HEADER

//------------------------------------------------
// max of vector
// max
// DEFINED IN HEADER

//------------------------------------------------
// square function
// sq
// DEFINED IN HEADER

//------------------------------------------------
// analogue of R function seq() for integers
vector<int> seq_int(int from, int to, int by) {
  int n = floor((to-from)/double(by)) + 1;
  vector<int> ret(n,from);
  for (int i=1; i<n; i++) {
    from += by;
    ret[i] = from;
  }
  return ret;
}

//------------------------------------------------
// Euclidian distance between points in 2 dimensions
// dist_euclid_2d
// DEFINED IN HEADER

//------------------------------------------------
// push back multiple values to vector
// push_back_multiple
// DEFINED IN HEADER

//------------------------------------------------
// erase particular element of a vector using efficient method. Warning - does
// not preserve original order of vector
// quick_erase
// DEFINED IN HEADER

//------------------------------------------------
// erase-remove idiom, for erasing all instances of a particular value from
// container
// erase_remove
// DEFINED IN HEADER

//------------------------------------------------
// test whether value can be found in vector
// is_in_vector
// DEFINED IN HEADER

//------------------------------------------------
// return unique values in a vector
// unique
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
// helper function for printing contents of a vector
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
  for (int i=0; i<n; i++) {
    Rcpp::Rcout << "*";
  }
  Rcpp::Rcout << "\n";
#else
  std::cout << title << " ";
  for (int i=0; i<n; i++) {
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
  return Rcpp::as<vector<bool> >(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<int> format.
vector<int> rcpp_to_vector_int(SEXP x) {
  return Rcpp::as<vector<int> >(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<double> format.
vector<double> rcpp_to_vector_double(SEXP x) {
  return Rcpp::as<vector<double> >(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<string> format.
vector<string> rcpp_to_vector_string(SEXP x) {
  return Rcpp::as<vector<string> >(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<bool>> format.
vector<vector<bool>> rcpp_to_matrix_bool(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<bool> > x_mat(nrow);
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<bool> >(x[i]);
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
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<int> >(x[i]);
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
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<double> >(x[i]);
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
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<string> >(x[i]);
  }
  return x_mat;
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
vector<vector<vector<int>>> rcpp_to_array_int(Rcpp::List x) {
  int n1 = int(x.size());
  vector< vector< vector<int> > > ret(n1);
  for (int i=0; i<n1; i++) {
    Rcpp::List x_i = x[i];
    int n2 = int(x_i.size());
    ret[i] = vector< vector<int> >(n2);
    for (int j=0; j<n2; j++) {
      ret[i][j] = Rcpp::as<vector<int> >(x_i[j]);
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
  vector< vector< vector<double> > > ret(n1);
  for (int i=0; i<n1; i++) {
    Rcpp::List x_i = x[i];
    int n2 = int(x_i.size());
    ret[i] = vector< vector<double> >(n2);
    for (int j=0; j<n2; j++) {
      ret[i][j] = Rcpp::as<vector<double> >(x_i[j]);
    }
  }
  return ret;
}
#endif

//------------------------------------------------
// read values from comma-separated text file to vector<int>
vector<int> file_to_vector_int(string file_path) {
  
  // initialise return object
  vector<int> ret;
  
  // read in values from comma-separated file
  ifstream infile(file_path);
  std::string line1, line2;
  int x;
  while (getline(infile, line1)) {
    istringstream ss(line1);
    while (getline(ss, line2, ',')) {
      if (line2.size() > 0) {
        istringstream(line2) >> x;
        ret.push_back(x);
      }
    }
  }
  
  return ret;
}

//------------------------------------------------
// read values from comma-separated text file to vector<double>
vector<double> file_to_vector_double(string file_path) {
  
  // initialise return object
  vector<double> ret;
  
  // read in values from comma-separated file
  ifstream infile(file_path);
  std::string line1, line2;
  double x;
  while (getline(infile, line1)) {
    istringstream ss(line1);
    while (getline(ss, line2, ',')) {
      if (line2.size() > 0) {
        istringstream(line2) >> x;
        ret.push_back(x);
      }
    }
  }
  
  return ret;
}

//------------------------------------------------
// read values from text file to vector<vector<double>>. Text file should be
// delimited by first line break, then comma. Lines do not all need to be same
// length, i.e. jagged matrices are allowed.
vector<vector<double>> file_to_matrix_double(string file_path) {
  
  // initialise return object
  vector<vector<double>> ret;
  
  // read in values from comma-separated file
  ifstream infile(file_path);
  std::string line1, line2, line3;
  double x;
  vector<double> v;
  while (getline(infile, line1)) {
    v.clear();
    istringstream ss(line1);
    while (getline(ss, line2, ',')) {
      if (line2.size() > 0) {
        istringstream(line2) >> x;
        v.push_back(x);
      }
    }
    ret.push_back(v);
  }
  
  return ret;
}

//------------------------------------------------
// calculate Cholesky decomposition of positive definite matrix sigma
void cholesky(vector<vector<double>> &chol, vector<vector<double>> &sigma) {
  
  for (int i = 0; i < int(sigma.size()); ++i) {
    for (int j = 0; j < (i+1); ++j) {
      chol[i][j] = sigma[i][j];
      if (i == j) {
        if (i > 0) {
          for (int k = 0; k < i; ++k) {
            chol[i][i] -= chol[i][k]*chol[i][k];
          }
        }
        chol[i][i] = sqrt(chol[i][i]);
      } else {
        if (j > 0) {
          for (int k = 0; k < j; ++k) {
            chol[i][j] -= chol[i][k]*chol[j][k];
          }
        }
        chol[i][j] /= chol[j][j];
      }
    }
  }
  
}
