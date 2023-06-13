// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// mcmc.cpp
list mcmc(const int chain, const bool burnin, const int iterations, const bool silent, const doubles_matrix<> theta_init, const strings theta_names, const integers transform_type, const doubles theta_min, const doubles theta_max, const integers infer_parameter, const list data, function ll_f, function lp_f, writable::list misc, const doubles_matrix<> proposal_sd_init, const integers_matrix<> acceptance_init, const double target_acceptance, const int swap, const doubles beta_init, const integers swap_acceptance_init, list blocks_list, const int n_unique_blocks, const int iteration_counter_init, cpp11::sexp rng_ptr);
extern "C" SEXP _drjacoby_mcmc(SEXP chain, SEXP burnin, SEXP iterations, SEXP silent, SEXP theta_init, SEXP theta_names, SEXP transform_type, SEXP theta_min, SEXP theta_max, SEXP infer_parameter, SEXP data, SEXP ll_f, SEXP lp_f, SEXP misc, SEXP proposal_sd_init, SEXP acceptance_init, SEXP target_acceptance, SEXP swap, SEXP beta_init, SEXP swap_acceptance_init, SEXP blocks_list, SEXP n_unique_blocks, SEXP iteration_counter_init, SEXP rng_ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(mcmc(cpp11::as_cpp<cpp11::decay_t<const int>>(chain), cpp11::as_cpp<cpp11::decay_t<const bool>>(burnin), cpp11::as_cpp<cpp11::decay_t<const int>>(iterations), cpp11::as_cpp<cpp11::decay_t<const bool>>(silent), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(theta_init), cpp11::as_cpp<cpp11::decay_t<const strings>>(theta_names), cpp11::as_cpp<cpp11::decay_t<const integers>>(transform_type), cpp11::as_cpp<cpp11::decay_t<const doubles>>(theta_min), cpp11::as_cpp<cpp11::decay_t<const doubles>>(theta_max), cpp11::as_cpp<cpp11::decay_t<const integers>>(infer_parameter), cpp11::as_cpp<cpp11::decay_t<const list>>(data), cpp11::as_cpp<cpp11::decay_t<function>>(ll_f), cpp11::as_cpp<cpp11::decay_t<function>>(lp_f), cpp11::as_cpp<cpp11::decay_t<writable::list>>(misc), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(proposal_sd_init), cpp11::as_cpp<cpp11::decay_t<const integers_matrix<>>>(acceptance_init), cpp11::as_cpp<cpp11::decay_t<const double>>(target_acceptance), cpp11::as_cpp<cpp11::decay_t<const int>>(swap), cpp11::as_cpp<cpp11::decay_t<const doubles>>(beta_init), cpp11::as_cpp<cpp11::decay_t<const integers>>(swap_acceptance_init), cpp11::as_cpp<cpp11::decay_t<list>>(blocks_list), cpp11::as_cpp<cpp11::decay_t<const int>>(n_unique_blocks), cpp11::as_cpp<cpp11::decay_t<const int>>(iteration_counter_init), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(rng_ptr)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_drjacoby_mcmc", (DL_FUNC) &_drjacoby_mcmc, 24},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_drjacoby(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
