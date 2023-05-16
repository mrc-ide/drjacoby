#include <cpp11.hpp>
#include <dust/r/random.hpp>

[[cpp11::register]]
double rngt(cpp11::sexp ptr, int chain){
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(ptr);
  auto& state = rng->state(chain - 1);
  double n1 = dust::random::random_real<double>(state);
  return n1;
}
