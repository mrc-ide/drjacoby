#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <cpp11.hpp>

double phi_to_theta(double phi, int transformation_type, double theta_min, double theta_max);
double theta_to_phi(double theta, int transformation_type, double theta_min, double theta_max);
double get_adjustment(double theta, double theta_prop, int transformation_type, double theta_min, double theta_max);

#endif
