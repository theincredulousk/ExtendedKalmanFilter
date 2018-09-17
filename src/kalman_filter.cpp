#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

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
  x_ = (F_ * x_);
  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = (z - H_ * x_);
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  // new state
  long dim = x_.size();
  x_ = (x_ + (K * y));
  P_ = (MatrixXd::Identity(dim, dim) - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  // cartesian back to polar
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float rho = sqrt(pow(px,2) + pow(py,2));
  float theta = atan2(py, px);

  // if rho was zero this goes to infinity, just picking a small value seems arbitrary...
  float rhod = (px * vx + py * vy);
  if(rho < 0.00001)
  {
    rho = 0.00001;
  }
  
  rhod = rhod / rho;

  VectorXd predicted = VectorXd(3);
  predicted << rho, theta, rhod;
  VectorXd delta = z - predicted;
  
  const float PI = M_PI;
  while(delta(1) > PI)
  {
    delta(1) -= 2.0*PI;
  }

  while(delta(1) < -1.0*PI)
  {
    delta(1) += 2.0*PI;
  }

  //std::cout << "delta " << delta << std::endl;

  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  long dim = x_.size();
  x_ = (x_ + (K * delta));
  P_ = (MatrixXd::Identity(dim, dim) - K * H_) * P_;

}
//RMSE: 0.106164,0.140782,0.433384,0.438315,
//AVG RMSE: 0.102325,0.138721,0.679333,0.532882,
