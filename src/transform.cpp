#include <cpp11.hpp>
#include "Rmath.h"

double phi_to_theta(double phi, int transformation_type, double theta_min, double theta_max) {

  double theta = phi; //default if transform type = 0
  if(transformation_type == 1){
    theta = theta_max - exp(phi);
  }
  if(transformation_type == 2){
    theta = exp(phi) + theta_min;
  }
  if(transformation_type == 3){
    theta = (theta_max * exp(phi) + theta_min) / (1 + exp(phi));
  }
  return(theta);
}

double theta_to_phi(double theta, int transformation_type, double theta_min, double theta_max) {

  double phi = theta; //default if transform type = 0
  if(transformation_type == 1){
    phi = log(theta_max - theta);
  }
  if(transformation_type == 2){
    phi = log(theta - theta_min);
  }
  if(transformation_type == 3){
    phi = log(theta - theta_min) - log(theta_max - theta);
  }

  return(phi);
}

double get_adjustment(double theta, double theta_prop, int transformation_type, double theta_min, double theta_max) {

  double adjustment;

  if(transformation_type == 0){
    adjustment = 0.0;
  }
  if(transformation_type == 1){
    adjustment = log(theta_max - theta_prop) - log(theta_max - theta);
  }
  if(transformation_type == 2){
    adjustment = log(theta_prop - theta_min) - log(theta - theta_min);
  }
  if(transformation_type == 3){
    adjustment = log(theta_max - theta_prop) + log(theta_prop - theta_min) - log(theta_max - theta) - log(theta - theta_min);

  }
  return(adjustment);
}
