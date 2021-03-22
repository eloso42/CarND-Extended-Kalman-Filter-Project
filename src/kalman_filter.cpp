#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

constexpr float PI = 3.14159265359f;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() :
    P_(MatrixXd::Identity(4,4)),
    F_(MatrixXd::Identity(4,4)),
    Q_(MatrixXd::Zero(4,4)),
    I_(MatrixXd::Identity(4,4)) {
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * predict the state
   */
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  // KF Measurement update step
  VectorXd y = z - H_ * x_;
  MatrixXd PHT = P_ * H_.transpose();
  MatrixXd S = H_ * PHT + R_;
  MatrixXd K = PHT * S.inverse();
  // new state
  x_ = x_ + K * y;
  P_ = (I_ - K*H_)*P_.transpose();
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  // KF Measurement update step
  VectorXd hx(3);
  hx(0) = sqrt(x_[0]*x_[0]+x_[1]*x_[1]);
  hx(1) = atan2(x_[1], x_[0]);
  hx(2) = (x_[0]*x_[2] + x_[1]*x_[3]) / hx(0);
  VectorXd y = z - hx;
  if (y(1) > PI) {
    y(1) -= 2*PI;
  }  else if (y(1) < -PI) {
    y(1) += 2*PI;
  }
  MatrixXd PHT = P_ * H_.transpose();
  MatrixXd S = H_ * PHT + R_;
  MatrixXd K = PHT * S.inverse();
  // new state
  x_ = x_ + K * y;
  P_ = (I_ - K*H_)*P_.transpose();
}
