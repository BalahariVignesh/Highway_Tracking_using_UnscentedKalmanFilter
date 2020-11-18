#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  
   std_a_ = 0.85;
   std_yawdd_ = 0.33;

  //The state dimension is 5 and augmented state dimension is 7
  n_x_=5;
  n_aug_ = 7;
  // create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  Xsig_pred_.fill(0.0);
  //define spreading parameter
  lambda_ = 3-n_aug_;//_aug or n_x
  
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.0);
  x_ = VectorXd(n_x_);
  P_ = MatrixXd(n_x_,n_x_);
  P_.fill(0.0);
  P_(0,0)=1;
  P_(1,1)=1;
  P_(2,2)=1;
  P_(3,3)=std_laspx_*std_laspx_;
  P_(4,4) =std_laspy_*std_laspy_;
  is_initialized_=false;
  time_us_=0.0;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  std::cout << "ProcessMeasurement Start" <<std::endl;
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_){
    time_us_ = meas_package.timestamp_;
    if(meas_package.sensor_type_==MeasurementPackage::RADAR)
    {
      // RADAR code
      double ro = meas_package.raw_measurements_[0];     
      double phi = meas_package.raw_measurements_[1];
      double ro_dot = meas_package.raw_measurements_[2];
      double vx = ro_dot*cos(phi);
      double vy = ro_dot*sin(phi);
      
      x_ << ro * cos(phi), ro * sin(phi), ro_dot, 0, 0;
       is_initialized_=true;
    }
    if (meas_package.sensor_type_==MeasurementPackage::LASER)
    {
      // LIDAR Code
      
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      is_initialized_=true;
    }
  
    return;
  }
  else{
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    
    Prediction(dt);
    if(meas_package.sensor_type_==MeasurementPackage::RADAR && use_radar_){
      UpdateRadar(meas_package);
    }
    if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
      UpdateLidar(meas_package);
    }
    time_us_ = meas_package.timestamp_;
  }
   std::cout << "ProcessMeasurement End " <<std::endl;
}

void UKF::Prediction(double delta_t) {

  std::cout << "Prediction Start" <<std::endl;
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  //NOTES for SELF
  // Generate Sigma points - Calculate Xsig
  // Augmentation - Calculate Xsig_aug
  // Sigma points Prediction - Calculate Xsig_pred
  // Predict mean and covariance - Calculate x(predicted state) and P(Predicted Covariance matrix)
  
  //Begin from augmentation part- Using x_ and P_ and generating Xsig_aug_
  // create augmented mean vector
  VectorXd x_aug_ = VectorXd(7);

  // create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(7, 7);

  // create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean state
    x_aug_.head(5)=x_;
    x_aug_(5)=0;
    x_aug_(6)=0;
    
  // create augmented covariance matrix
    P_aug_.fill(0);
    P_aug_.topLeftCorner(5,5)=P_;
    P_aug_(5,5)=std_a_*std_a_;
    P_aug_(6,6)=std_yawdd_*std_yawdd_;
  // create square root matrix
    MatrixXd L = P_aug_.llt().matrixL();
    
  // create augmented sigma points
    Xsig_aug_.col(0) = x_aug_;
    for(int i = 0;i<n_aug_;i++){
        Xsig_aug_.col(i+1)=x_aug_+sqrt(lambda_+n_aug_)*L.col(i);
        Xsig_aug_.col(i+1+n_aug_)=x_aug_-sqrt(lambda_+n_aug_)*L.col(i);
    }
  
  //Begin Sigma Point prediction part - Using X_sig_aug_ and generating Xsig_pred_
  for(int i=0;i<2*n_aug_+1;++i){
        double p_x = Xsig_aug_(0,i);
        double p_y = Xsig_aug_(1,i);
        double v = Xsig_aug_(2,i);
        double yaw = Xsig_aug_(3,i);
        double yawd = Xsig_aug_(4,i);
        double nu_a = Xsig_aug_(5,i);
        double nu_yawdd = Xsig_aug_(6,i);
        
        double px_p, py_p;
        
        //avoid division by zero
        if(fabs(yawd)>0.001){
            px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t)-sin(yaw));
            py_p = p_y + v/yawd * (cos(yaw)-cos(yaw+yawd*delta_t));
        } else{
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        //write predicted sigma point into right column
        
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }

    //Begin predict mean and covariance - Using Xsig_pred_ and generating x and P
    // create vector for weights
      //VectorXd weights = VectorXd(2*n_aug_+1);
      
      

      // // set weights
        double weight_0 = lambda_/(lambda_+n_aug_);
        weights_(0) = weight_0;
        //++i changed to i++
        for(int i =1; i<2*n_aug_+1;i++){
            double weight = 0.5/(n_aug_+lambda_);
            weights_(i)=weight;
        }
        
      // predict state mean
        x_.fill(0.0);
        for (int i = 0;i<2*n_aug_ + 1;i++){
            x_=x_+weights_(i)*Xsig_pred_.col(i);
        }
      // predict state covariance matrix
        P_.fill(0.0);
        for(int i = 0; i<2*n_aug_ + 1;i++){
            VectorXd x_diff = Xsig_pred_.col(i)-x_;
            //normalizing angles.
            while(x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
            while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
            
            P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
        }
      
      // print result
      // std::cout << "Predicted state" << std::endl;
      // std::cout << x << std::endl;
      // std::cout << "Predicted covariance matrix" << std::endl;
      // std::cout << P << std::endl;
      std::cout << "Prediction End" <<std::endl;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {

  std::cout << "UpdateLidar Start" <<std::endl;
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  //Predict Lidar Measurements
  // set measurement dimension, lidar has only 2 dimensions unlike radar
  int n_z = 2;
  // create matrix for sigma points in measurement space 
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);//Check declaration
    Zsig.fill(0.0);
    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    
    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
  // transform sigma points into measurement space
      for(int i = 0; i<Xsig_pred_.cols();i++){
         
          Zsig.col(i) << Xsig_pred_(0,i),Xsig_pred_(1,i);

      }
    
    // calculate mean predicted measurement
    z_pred.fill(0.0);
    for (int i =0;i<Zsig.cols();i++){

        z_pred = z_pred + weights_(i)*Zsig.col(i);
    }
    // calculate innovation covariance matrix S
    S.fill(0.0);
    for (int i=0;i<Zsig.cols();i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        S = S + weights_(i) * z_diff * z_diff.transpose();
        }
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R<<std_laspx_*std_laspx_,0,
      0,std_laspy_*std_laspy_;
    S = S+R;
  //UKF Update
  // create matrix for cross correlation Tc
      MatrixXd Tc = MatrixXd(n_x_, n_z);
    //Incoming radar measurement assigned to z
    VectorXd z = VectorXd(n_z);
    z<<meas_package.raw_measurements_(0),meas_package.raw_measurements_(1);
    // calculate cross correlation matrix
    Tc.fill(0.0);
    for(int i = 0; i<Zsig.cols();i++){
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while(z_diff(1)< -M_PI) z_diff(1)+=2.*M_PI;
        
        //state difference
        VectorXd x_diff = Xsig_pred_.col(i)-x_;
        //angle normalization
        while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
        while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    // calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();
  std::cout << "2  " <<std::endl;
    //residual
    VectorXd z_diff = z-z_pred;

    // //angle normalization
    while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)< -M_PI) z_diff(1)=+2.*M_PI;

    // update state mean and covariance matrix
    x_= x_+K*z_diff;
    P_=P_-K*S*K.transpose();
  std::cout << " 3 " <<std::endl;
    //optional
    //double nis_lidar = z_diff.transpose() * S.inverse() * z_diff;
    std::cout << "UpdateLidar  End" <<std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

  std::cout << "UpdateRadar Start" <<std::endl;
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  // NOTES TO SELF
  // Predict Radar Measurements 
  //   Transform the sigma points in measurement space
  //   calculate mean predicted Measurement
  //   calculate innovation covariance matrix S
  //   add measurement noise covariance matrix
  // UKF update


  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  // create matrix for sigma points in measurement space 
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);//Check declaration

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    
    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);


    // transform sigma points into measurement space
      for(int i = 0; i<2*n_aug_+1;++i){
          double p_x = Xsig_pred_(0,i);
          double p_y = Xsig_pred_(1,i);
          double v = Xsig_pred_(2,i);
          double yaw = Xsig_pred_(3,i);
          
          double v1 = cos(yaw)*v;
          double v2 = sin(yaw)*v;
          double rho_d;
          //condition to division b zero
          double rho = sqrt(p_x*p_x+p_y*p_y);
          if(rho <0.001){
            rho_d = (p_x*v1 + p_y*v2) / 0.001;
          }else{
            rho_d = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x+p_y*p_y);
          }
          
          Zsig(0,i) = sqrt(p_x*p_x+p_y*p_y);
          Zsig(1,i) = atan2(p_y,p_x);
          Zsig(2,i) = rho_d;
      }
    
    // calculate mean predicted measurement
    z_pred.fill(0.0);
    for (int i =0;i<2*n_aug_+1;++i){
        z_pred = z_pred + weights_(i)*Zsig.col(i);
    }
    
    // calculate innovation covariance matrix S
    S.fill(0.0);
    for (int i=0;i<2*n_aug_+1;i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        // angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      S = S + weights_(i) * z_diff * z_diff.transpose();
    
    }
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R<<std_radr_*std_radr_,0,0,
      0,std_radphi_*std_radphi_,0,
      0,0,std_radrd_*std_radrd_;
    S = S+R;

   

    // print result
    // std::cout << "z_pred: " << std::endl << z_pred << std::endl;
    // std::cout << "S: " << std::endl << S << std::endl;

    //UKF Update
    // create matrix for cross correlation Tc
      MatrixXd Tc = MatrixXd(n_x_, n_z);
    //Incoming radar measurement assigned to z
    VectorXd z = VectorXd(n_z);
    z<<meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),meas_package.raw_measurements_(2);
      // calculate cross correlation matrix
        Tc.fill(0.0);
        for(int i = 0; i<2*n_aug_ + 1;i++){
            //residual
            VectorXd z_diff = Zsig.col(i) - z_pred;
            //angle normalization
            while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
            while(z_diff(1)< -M_PI) z_diff(1)+=2.*M_PI;
            
            //state difference
            VectorXd x_diff = Xsig_pred_.col(i)-x_;
            //angle normalization
            while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
            while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
            
            Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
        }
      // calculate Kalman gain K;
        MatrixXd K = Tc * S.inverse();
        
        //residual
        VectorXd z_diff = z-z_pred;
        
        //angle normalization
        while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while(z_diff(1)< -M_PI) z_diff(1)=+2.*M_PI;
        
      // update state mean and covariance matrix
        x_= x_+K*z_diff;
        P_=P_-K*S*K.transpose();
      

      // print result
      // std::cout << "Updated state x: " << std::endl << x << std::endl;
      // std::cout << "Updated state covariance P: " << std::endl << P << std::endl;
    //optional
    //double nis_radar = z_diff.transpose() * S.inverse() * z_diff;
  std::cout << "UpdateRadar End" <<std::endl;
}