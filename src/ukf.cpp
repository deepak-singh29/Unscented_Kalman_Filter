#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = x_.size();
  n_aug_ = n_x_ + 2;

  lambda_ = 3 - n_x_;
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0)=lambda_/(lambda_+n_aug_);
  for(int i =1;i<weights_.size();i++)
  {
    weights_(i)=1/(2*(lambda_+n_aug_));
  }
  
  is_initialized_ = false;
  time_us_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(!is_initialized_){
	//cout<<"UKF Initialisation"<<endl;
	//x_ << 1,1,1,1,1;
	  
	P_ = MatrixXd(5, 5);
	P_ << 1,0,0,0,0,
		  0,1,0,0,0,
		  0,0,1,0,0,
		  0,0,0,1,0,
		  0,0,0,0,1;
	
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho,phi,rho_dot,px,py,v;
      rho = meas_package.raw_measurements_(0);
      phi = meas_package.raw_measurements_(1);
      rho_dot = meas_package.raw_measurements_(2);
      px = rho * cos(phi);
      py = rho * sin(phi);
      v = sqrt(rho_dot*cos(phi) * rho_dot*cos(phi) + rho_dot*sin(phi) * rho_dot*sin(phi));

      x_ << px,py,v,0.0,0.0;
      
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),0.0,0.0,0.0;
      
    }
	time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  //cout << x_<<endl;
  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  //cout <<"time difference :"<<dt<<endl;
  Prediction(dt);
  //cout<<"predict call"<<endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_)) {
    // Radar updates
    // ekf_.R_ = R_radar_;
    UpdateRadar(meas_package);
	//cout<<"update Radar call"<<endl;
  } else  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_){
    // Laser updates
    // ekf_.R_ = R_laser_;
	//cout<<"update Radar before call"<<endl;
    UpdateLidar(meas_package);
	//cout<<"update Radar call"<<endl;
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  //--------------------------------------------------
  // Sigma Point Generation
  //--------------------------------------------------
  
  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  
  // generating Sigma points
  MatrixXd SigH1,SigH2 = MatrixXd(n_x_, n_x_);
  MatrixXd t = MatrixXd(1,5);
  t.setOnes();
  MatrixXd addS = x_ * t;
  SigH1 = addS + sqrt(n_x_ + lambda_) * A;
  SigH2 = addS - sqrt(n_x_ + lambda_) * A;
  Xsig <<  x_,SigH1,SigH2;
  //cout<< Xsig.size() <<endl;
  
  //--------------------------------------------------
  // Sigma Point Augmentation
  //--------------------------------------------------
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state
  x_aug << x_,
            0,
            0;
  //create augmented covariance matrix
    P_aug.setZero();
    P_aug.topLeftCorner(5, 5) = P_;
    P_aug.bottomRightCorner(2, 2) << std_a_ * std_a_,0,
                                    0, std_yawdd_ * std_yawdd_ ;
  //create square root matrix
  MatrixXd AAug = P_aug.llt().matrixL();
  
  //create augmented sigma points
  
  MatrixXd SigH1_aug,SigH2_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd t_a = MatrixXd(1,n_aug_);
  t_a.setOnes();
  MatrixXd addS_aug = x_aug * t_a;
  SigH1_aug = addS_aug + sqrt(n_aug_ + lambda_) * AAug;
  SigH2_aug = addS_aug - sqrt(n_aug_ + lambda_) * AAug;
  
  Xsig_aug << x_aug,SigH1_aug,SigH2_aug;
  //--------------------------------------------------
  
  
  //--------------------------------------------------
  // Sigma Point prediction
  //--------------------------------------------------
  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  for(int i = 0;i<Xsig_aug.cols();i++)
  {
    
    double psi_dot_k = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_psi = Xsig_aug(6,i);
    double v_k = Xsig_aug(2,i);
    double psi_k = Xsig_aug(3,i);
    //   Position vector x and y component
    if(psi_dot_k != 0)
    {
        Xsig_pred_(0,i) =  Xsig_aug(0,i) + ((v_k/psi_dot_k)*( sin(psi_k + psi_dot_k * delta_t) - sin(psi_k))) + (0.5 * delta_t * delta_t * cos(psi_k) * nu_a);
        Xsig_pred_(1,i) =  Xsig_aug(1,i) + ((v_k/psi_dot_k)*( -cos(psi_k + psi_dot_k * delta_t) + cos(psi_k))) + (0.5 * delta_t * delta_t * sin(psi_k) * nu_a);
    }
    else
    {
        Xsig_pred_(0,i) =  Xsig_aug(0,i) + ((v_k * cos(psi_k) * delta_t) + (0.5 * delta_t * delta_t * cos(psi_k) * nu_a));
        Xsig_pred_(1,i) =  Xsig_aug(1,i) + ((v_k * sin(psi_k) * delta_t) + (0.5 * delta_t * delta_t * sin(psi_k) * nu_a));
    }
    //   Velocity component
    Xsig_pred_(2,i) = v_k + delta_t * nu_a;
    //   Yaw angle component
    Xsig_pred_(3,i) = psi_k + (psi_dot_k * delta_t) + (0.5 * delta_t * delta_t * nu_psi);
    //   Yaw angle change rate component
    Xsig_pred_(4,i) = psi_dot_k + (delta_t * nu_psi);
  }
  //--------------------------------------------------
  
  //--------------------------------------------------
  // Predicted Mean Calculation
  //--------------------------------------------------
  //predict state mean
  x_.fill(0.0);
  for(int j=0;j<weights_.size();j++)
  {
    x_ = x_ + weights_(j) * Xsig_pred_.col(j);
  }
  //predict state covariance matrix
   P_.fill(0.0);
   for(int i=0;i<weights_.size();i++)
  {
    MatrixXd x_dif = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_dif(3)> M_PI) x_dif(3)-=2.*M_PI;
    while (x_dif(3)<-M_PI) x_dif(3)+=2.*M_PI;
    P_ = P_ + weights_(i) * x_dif * x_dif.transpose();
  }
  //--------------------------------------------------
  //cout<<"X_Pred :"<<x_<<endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //cout << "meas_package.raw_measurements_" << meas_package.raw_measurements_<<endl;
  VectorXd z = meas_package.raw_measurements_;
  MatrixXd H_ = MatrixXd(2,5);
  H_ << 1,0,0,0,0,
		0,1,0,0,0;
  MatrixXd R_ = MatrixXd(2,2);
  R_ << std_laspx_* std_laspx_,0,
		0,std_laspy_*std_laspy_;
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
	NIS_laser = y.transpose() * Si * y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //transform sigma points into measurement space
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
   //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_radar = z_diff.transpose() * S.inverse() * z_diff;
  cout<<NIS_radar<<","<<endl;
}
