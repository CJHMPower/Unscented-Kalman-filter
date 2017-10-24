#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ =1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_=5;
  n_aug_=7;
  lambda_ =3-n_aug_;
  is_initialized_ = false;

  //initial time step
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
  if(is_initialized_==false){
    double px;
    double py;
    double v;
    double phi;
    double phi_dot;
    weights_=VectorXd(2*n_aug_+1);
    weights_(0)=float(lambda_)/(lambda_+n_aug_);
    for(int i=1;i<2*n_aug_+1;i++){
      weights_(i)=1.0/(2*(lambda_+n_aug_));
    }
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    double r=meas_package.raw_measurements_(0);
    double theta=meas_package.raw_measurements_(1);
    double r_dot=meas_package.raw_measurements_(2);
    px=r*cos(theta);
    py=r*sin(theta);
    v=fabs(r_dot);
    phi=0;
    phi_dot=0;
    }
    else if(meas_package.sensor_type_==MeasurementPackage::LASER){
    px=meas_package.raw_measurements_(0);
    py=meas_package.raw_measurements_(1);
    
    double r=sqrt(px*px+py*py);
    
    v=0;
    phi=0;
    phi_dot=0;
    }
    is_initialized_=true;
    x_=VectorXd(n_x_);
    x_<<px,py,v,phi,phi_dot;
   
    P_=MatrixXd::Identity(n_x_, n_x_);
    time_us_=meas_package.timestamp_;
    return;
  }
  double dt=(meas_package.timestamp_-time_us_)/1000000.0;

  time_us_=meas_package.timestamp_;
  
  Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::LASER){
     UpdateLidar(meas_package);
  }
  else if(meas_package.sensor_type_==MeasurementPackage::RADAR){
     UpdateRadar(meas_package);
  }
  return;
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
  Xsig_pred_=MatrixXd(n_x_, 2*n_aug_+1);
  MatrixXd Xsig_aug_=MatrixXd(n_aug_, 2*n_aug_+1);
  VectorXd x_aug=VectorXd(n_aug_);
  MatrixXd P_aug=MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  x_aug.fill(0.0);
  x_aug.head(n_x_)=x_;
  x_aug(n_x_)=0;
  x_aug(n_x_+1)=0;
  
  P_aug.topLeftCorner(n_x_,n_x_)=P_;
  P_aug(n_x_,n_x_)=std_a_*std_a_;
  P_aug(n_x_+1, n_x_+1)=std_yawdd_*std_yawdd_;
  Xsig_aug_.col(0)=x_aug;
  
  MatrixXd A=P_aug.llt().matrixL();
  for(int i=0;i<n_aug_;i++){
      Xsig_aug_.col(i+1)=x_aug+sqrt(lambda_+n_aug_)*A.col(i);
      Xsig_aug_.col(i+1+n_aug_)=x_aug-sqrt(lambda_+n_aug_)*A.col(i);
  }
 
  
  double delta_t2=delta_t*delta_t;
  for(int i=0;i<2*n_aug_+1;i++){
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);
    // Precalculate sin and cos for optimization
    double sin_yaw = sin(yaw);
    double cos_yaw = cos(yaw);
    double arg = yaw + yawd*delta_t;
    
    // Predicted state values
    double px_p, py_p;
    // Avoid division by zero
    if (fabs(yawd) > 0.0001) {	
	double v_yawd = v/yawd;
        px_p = p_x + v_yawd * (sin(arg) - sin_yaw);
        py_p = p_y + v_yawd * (cos_yaw - cos(arg) );
    }
    else {
	double v_delta_t = v*delta_t;
        px_p = p_x + v_delta_t*cos_yaw;
        py_p = p_y + v_delta_t*sin_yaw;
    }
    double v_p = v;
    double yaw_p = arg;
    double yawd_p = yawd;

    // Add noise
    px_p += 0.5*nu_a*delta_t2*cos_yaw;
    py_p += 0.5*nu_a*delta_t2*sin_yaw;
    v_p += nu_a*delta_t;
    yaw_p += 0.5*nu_yawdd*delta_t2;
    yawd_p += nu_yawdd*delta_t;

    // Write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  

  
  x_ = Xsig_pred_ * weights_;
  P_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
      MatrixXd cov=MatrixXd(n_x_,n_x_);
      VectorXd diff=VectorXd(n_x_);
      diff=Xsig_pred_.col(i)-x_; 
      while(diff(3)<-M_PI) diff(3)+=2*M_PI;
      while(diff(3)>M_PI) diff(3)-=2*M_PI;
      cov=diff*diff.transpose();
      P_=P_+weights_(i)*cov;
  } 
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
  int n_z=2;
  MatrixXd z_sig=MatrixXd(n_z,2*n_aug_+1);
  VectorXd z_pred=VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z); 
  z_pred.fill(0.0);
  S.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
    z_sig(0,i)=Xsig_pred_(0,i);
    z_sig(1,i)=Xsig_pred_(1,i);
    z_pred+=weights_(i)*z_sig.col(i);
  }
  MatrixXd R=MatrixXd(n_z,n_z);
  R<<std_laspx_*std_laspx_,0,
  0, std_laspy_*std_laspy_;
  
  for(int i=0;i<2*n_aug_+1;i++){
      VectorXd diff=z_sig.col(i)-z_pred;
      MatrixXd cov=weights_(i)*(diff*diff.transpose());
      S+=cov;
  }
  S=S+R;
  
  VectorXd z=meas_package.raw_measurements_;
  MatrixXd Tc = MatrixXd(n_x_, n_z);  
  Tc.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
      VectorXd diff_x=Xsig_pred_.col(i)-x_;
      while(diff_x(3)>M_PI) diff_x(3)-=2*M_PI;
      while(diff_x(3)<-M_PI) diff_x(3)+=2*M_PI;
      VectorXd diff_z=z_sig.col(i)-z_pred;
      MatrixXd update=weights_(i)*(diff_x*diff_z.transpose());
      
      Tc+=update;
  }
  
 
  MatrixXd K=Tc*S.inverse();
  VectorXd z_diff=z-z_pred;
  x_=x_+K*(z-z_pred);
  P_=P_-K*S*K.transpose();
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  return;
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z=3;
  MatrixXd z_sig=MatrixXd(n_z,2*n_aug_+1);
  VectorXd z_pred=VectorXd(n_z);
  z_pred.fill(0.0);
  MatrixXd S = MatrixXd(n_z,n_z); 
  S.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
    double px=Xsig_pred_(0,i);
    double py=Xsig_pred_(1,i);
    double vk=Xsig_pred_(2,i);
    double phi=Xsig_pred_(3,i);
    double r=sqrt(px*px+py*py);
    z_sig(0,i)=r;
    double yaw=atan2(py,px);
    while(yaw>M_PI) yaw-=2.*M_PI;
    while(yaw<-M_PI) yaw+=2.*M_PI;
    z_sig(1,i)=yaw;
    if(fabs(r)>0.0001){
    z_sig(2,i)=(px*cos(phi)*vk+py*sin(phi)*vk)/r;
    }
    else{
    z_sig(2,i)=0;
    }
   
  }
  for(int i=0;i<2*n_aug_+1;i++){
     z_pred=z_pred+weights_(i)*z_sig.col(i);
  }
 
  MatrixXd R=MatrixXd(n_z,n_z);
  R<<std_radr_*std_radr_ , 0 , 0,
  0, std_radphi_*std_radphi_,0,
  0,0,std_radrd_*std_radrd_;  

  for(int i=0;i<2*n_aug_+1;i++){
      VectorXd diff=z_sig.col(i)-z_pred;
      while(diff(1)>M_PI) diff(1)-=2.*M_PI;
      while(diff(1)<-M_PI) diff(1)+=2.*M_PI;
      MatrixXd cov=weights_(i)*diff*diff.transpose();
      S+=cov;
  }
  S=S+R;

  VectorXd z=meas_package.raw_measurements_;
  MatrixXd Tc = MatrixXd(n_x_, n_z);  
  Tc.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
      VectorXd diff_x=Xsig_pred_.col(i)-x_;
      VectorXd diff_z=z_sig.col(i)-z_pred;
      while(diff_z(1)>M_PI) diff_z(1)-=2.*M_PI;
      while(diff_z(1)<-M_PI) diff_z(1)+=2.*M_PI;
      while(diff_x(3)>M_PI) diff_x(3)-=2.*M_PI;
      while(diff_x(3)<-M_PI) diff_x(3)+=2.*M_PI;
      Tc+=weights_(i)*diff_x*diff_z.transpose();
  }
  
 
  MatrixXd K=Tc*S.inverse();
  VectorXd diff_z=z-z_pred;
  x_=x_+K*(z-z_pred);
  
  P_=P_-K*S*K.transpose();
  NIS_radar_ = diff_z.transpose() * S.inverse() * diff_z;
  return;
}
