#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
    return rmse;
  }

  // accumulate squared residuals
  for (int i=0; i < estimations.size(); ++i) {
    // ... your code here
    VectorXd diff = estimations[i] - ground_truth[i];
    VectorXd sqdiff = diff.array() * diff.array();
    rmse += sqdiff;
  }

  // calculate the mean
  rmse /= estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  Hj.fill(0.0f);

  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // TODO: YOUR CODE HERE
  float pxy2 = px*px + py*py;
  float srpxy2 = sqrtf(pxy2);
  float pxy23 = pxy2 * srpxy2;
  if (fabs(pxy2) > 1e-6f) {
    Hj(0,0) = px / srpxy2;
    Hj(0,1) = py / srpxy2;
    Hj(1,0) = -py / pxy2;
    Hj(1,1) = px / pxy2;
    Hj(2,0) = py * (vx *py - vy * px) / pxy23;
    Hj(2,1) = px * (vy * px - vx * py) / pxy23;
    Hj(2,2) = Hj(0,0);
    Hj(2,3) = Hj(0,1);
  }
  return Hj;
}
