#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

    if(estimations.size() == 0 || estimations.size() != ground_truth.size())
    {
        cout << "error" << std::endl;
        return rmse;
    }

    VectorXd residual(4);
    residual << 0,0,0,0;
    
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        VectorXd sample = estimations[i] - ground_truth[i];
        //cout << "sample" << endl << sample << endl;
        sample = pow(sample.array(), 2);
        //cout << "sample2" << endl << sample << endl;
		residual += sample;
		//cout << "residual" << endl << residual << endl;
	}

	//calculate the mean
	//cout << "est size " << estimations.size() << endl;
	//cout << "est inv " << 1.0/estimations.size() << endl;
    residual = (1.0/estimations.size()) * residual.array();
    //cout << "residual mean" << endl << residual << endl;

	//calculate the squared root
	rmse = residual.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3,4);
	Hj << 	0,0,0,0,
			0,0,0,0,
			0,0,0,0;
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE 
	float sum_pxpy = pow(px,2) + pow(py,2);
    float sqrt_pxpy = sqrt(sum_pxpy);
	//check division by zero
	if(sum_pxpy == 0.0 || sqrt_pxpy == 0.0)
	{
	    cout << "err divide by 0" << endl;
	    //sum_pxpy = 0.0001;
	    //sqrt_pxpy = 0.01;
	    return Hj;
	}
	//compute the Jacobian matrix
    Hj(0,0) = px / sqrt_pxpy;
	Hj(0,1) = py / sqrt_pxpy;
	
	Hj(1,0) = -1 * (py / sum_pxpy);
	Hj(1,1) = 1 * (px / sum_pxpy);
	
	Hj(2,0) = (py * (vx*py - vy*px)) / pow(sum_pxpy, (1.5));
	Hj(2,1) = (px * (vy*px - vx*py)) / pow(sum_pxpy, (1.5));
	Hj(2,2) = px / sqrt_pxpy;
	Hj(2,3) = py / sqrt_pxpy;;
	return Hj;
}
