#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  noise_ax_ = 9.0;
  noise_ay_ = 9.0;
  H_laser_ << 1, 0, 0, 0,
             0, 1, 0, 0;

  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 10, 0, 0, 0,
              0, 100, 0, 0,
              0, 0, 500, 0,
              0, 0, 0, 1500;
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << 0,0,0,0,
             0,0,0,0,
             0,0,0,0,
             0,0,0,0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1.5, 0.6, 7, 0.60; // initialized velocity low - we have no delta t yet?
    //ekf_.x_ << 0.3, 0.3, 0, 0; 
    previous_timestamp_ = measurement_pack.timestamp_;
    VectorXd p_raw = measurement_pack.raw_measurements_;
    cout << "raw INIT = " << p_raw << endl;
    //cout << "EKF init " << p_raw;

    //float v_est = 

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = p_raw(0);
      float theta = p_raw(1);
      //float rho_dot = p_raw(2);

      ekf_.x_(0) = rho * cos(theta);
      ekf_.x_(1) = rho * sin(theta);

      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            ekf_.x_(0) = p_raw(0);
            ekf_.x_(1) = p_raw(1);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  cout << "raw INPUT = " << measurement_pack.raw_measurements_ << endl;
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //cout << "Timestamp= " << measurement_pack.timestamp_ << std::endl;
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  // Refer to Laser Measurements Part 3
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_  << (pow(dt, 4) / 4)*noise_ax_,    0, (pow(dt, 3) / 2)*noise_ax_ , 0
              ,0,(pow(dt, 4) / 4)*noise_ay_, 0, (pow(dt, 3) / 2)*noise_ay_
              ,(pow(dt, 3) / 2)*noise_ax_,   0, pow(dt, 2)*noise_ax_       , 0
              ,0, (pow(dt, 3)/2)*noise_ay_,  0, pow(dt, 2)*noise_ay_;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    //cout << "ekfH R" << ekf_.H_ << endl;
    ekf_.R_ = R_radar_;
    //cout << "ekfR R" << ekf_.R_ << endl;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
  else 
  {
    ekf_.H_ = H_laser_;
    //cout << "ekfH L" << ekf_.H_ << endl;
    ekf_.R_ = R_laser_;
    //cout << "ekfR L" << ekf_.R_ << endl;
    
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ output = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
