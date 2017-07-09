/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

// declare a random engine to be used across multiple and various method calls
//static default_random_engine gen;

//const int N_particles = 100;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    // initialize random engine
	default_random_engine gen;
	//initialize number of particles 
	num_particles = 100;

	//Set standard deviations for x, y, and theta.
	double std_x     = std[0]; 
	double std_y     = std[1]; 
	double std_theta = std[2];

	//define normal distribution for sensor noise (zero mean and standard deviation std)
	normal_distribution<double> noise_dist_x (0, std_x);
	normal_distribution<double> noise_dist_y (0, std_y);
	normal_distribution<double> noise_dist_theta (0, std_theta);

  	// initialize particles
	for (int i = 0; i < num_particles; i++) 
	{		
		Particle  p;
		p.id = i;
		p.x = x;
		p.y = y;
		p.theta = theta;
		p.weight = 1.0;

		//add noise 
		p.x 	+= noise_dist_x(gen);
		p.y 	+= noise_dist_y(gen);
		p.theta += noise_dist_theta(gen);
		

		particles.push_back(p);
		weights.push_back(p.weight);

		//print samples 
		cout << "initial position x in GPS:" << p.x << endl;
		cout << "initial position y in GPS:" << p.y << endl;
		cout << "initial orientation :" << p.theta << endl;
	}
  	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// initialize random engine
	default_random_engine gen;

	//Set standard deviations for x, y, and theta.
	double std_x     = std_pos[0]; 
	double std_y     = std_pos[1]; 
	double std_theta = std_pos[2];

	//define normal distribution for sensor noise (zero mean and standard deviation std)
	normal_distribution<double> noise_dist_x(0, std_x );
	normal_distribution<double> noise_dist_y(0, std_y );
	normal_distribution<double> noise_dist_theta(0, std_theta );

	for (int i = 0; i < num_particles; i++)
	{
		Particle &p = particles[i];

		if (fabs(yaw_rate) < 1e-5)
		{
			p.x += velocity * cos(p.theta) * delta_t;
			p.y += velocity * sin(p.theta) * delta_t;
		}
		else
		{
			p.x += ( velocity / yaw_rate ) * ( sin(p.theta + yaw_rate * delta_t) - sin(p.theta) );
			p.y += ( velocity / yaw_rate ) * ( cos(p.theta) - cos(p.theta + yaw_rate * delta_t) );
			p.theta += yaw_rate * delta_t; 
		}	
		//add noise 
		p.x 	+= noise_dist_x(gen);
		p.y 	+= noise_dist_y(gen);
		p.theta += noise_dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	/*
	struct LandmarkObs 
	{
	int id;				// Id of matching landmark in the map.
	double x;			// Local (vehicle coordinates) x position of landmark observation [m]
	double y;			// Local (vehicle coordinates) y position of landmark observation [m]
	};
	*/

	for (int  i = 0; i < observations.size(); i++)
	{
		// intialize minimum distance to maximum possible
		double min_distance = numeric_limits<double>::max();

		// current observation
		LandmarkObs &obs = observations[i];

		for (int j = 0; j < predicted.size(); j++)
		{
			// current prediction
			LandmarkObs pred = predicted[j];

			// compute the distance between current predicted and landmark
			double distance = dist(obs.x, obs.y, pred.x, pred.y);

			// find the predicted nearest to the current landmark
			if (distance < min_distance) 
			{
				min_distance = distance;

				// set the current observation to the nearest predicted (predicted's id)
				obs.id = j;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	/*
	class Map {
	public:
	
	struct single_landmark_s{
		int id_i ; // Landmark ID
		float x_f; // Landmark x-position in the map (global coordinates)
		float y_f; // Landmark y-position in the map (global coordinates)
	};

	std::vector<single_landmark_s> landmark_list ; // List of landmarks in the map

	};
	*/

	// for each particle
	for (int i = 0; i <  num_particles; i++)
	{
		Particle &p = particles[i];
		//get the particle x,y coordinates (global coordinates)
		double p_x 		= p.x;
		double p_y 	 	= p.y;
		double p_theta 	= p.theta;


		// create a copy of a list of observations from vehicle coordinates to map coodinates (global coordinates)
		vector<LandmarkObs> transformed_obs;

		//transform observation to map coordinates (global coordinates)
		for (int j = 0; j < observations.size(); j++)
		{
			LandmarkObs obs = observations[j];
			LandmarkObs transformed_obs_landmark;

			//  rotation AND translation transformation transform_x = T * x 
			// using homogeneous transformation matrix T = [ cos(theta), -sin(theta), x_t
			//												 sin(theta), cos(theta) , y_t
			//														 0, 		   0, 1 ]
			//
			// x = [observation.x, observation.y, observation.theta]
			// x_t = p_x
			// y_t = p_y

			// roate and translate
			transformed_obs_landmark.x = cos(p_theta) * obs.x - sin(p_theta) * obs.y + p_x;
			transformed_obs_landmark.y = sin(p_theta) * obs.x + cos(p_theta) * obs.y + p_y;
			//transformed_obs_landmark.id = obs.id;

			// add transform landmark to a vector
			transformed_obs.push_back(transformed_obs_landmark);
		}

		// predict measurements
		//create a vector for the map landmark locations predicted within the sensor range of the particle
		vector<LandmarkObs> predicted_landmarks;
		
		// for each map landmark
		for (int k = 0; k < map_landmarks.landmark_list.size(); k++){

			// get map landmark id, x, y coordinates  (global coordinates)
			double lm_x = map_landmarks.landmark_list[k].x_f;
			double lm_y = map_landmarks.landmark_list[k].y_f;
			int lm_id = map_landmarks.landmark_list[k].id_i;

			// find landmarks within the sensor range of the particle considering a circular region
			if (dist(lm_x, lm_y, p_x, p_y) <= sensor_range) 
			// find landmarks within the sensor range of the particle considering a rectangular region
			//if (fabs(lm_x - p_x) <= sensor_range && fabs(lm_y - p_y) <= sensor_range)
			{
				// create predicted landmark 
				LandmarkObs pred_landmark;

				pred_landmark.id = lm_id;
				pred_landmark.x = lm_x;
				pred_landmark.y = lm_y;

				// add predicted to vector predicted_landmarks
				predicted_landmarks.push_back(pred_landmark);
			}
		}

		// perform dataAssociation for predicted_landmarks and transformed_obs on current particle
		dataAssociation(predicted_landmarks, transformed_obs);

		// 	update particle weight
		// initialize the weight
		p.weight = 1.0;

		for ( int l = 0; l < transformed_obs.size(); l++){

			LandmarkObs obs = transformed_obs[l];
			LandmarkObs pred = predicted_landmarks[obs.id];

			// calculate the weight for this observation with multi-variant Gaussian 
			double std_x = std_landmark[0];
			double std_y = std_landmark[1];

			double obs_weight = (1/ (2 * M_PI * std_x * std_y) ) * exp(- ( pow(pred.x-obs.x,2) / (2 * pow(std_x,2)) + pow(pred.y-obs.y,2) / (2 * pow(std_y,2)) ) );

			//calculate the total observation weight
			p.weight *= obs_weight;				
		}
		weights[i] = p.weight;
	}
}


		/*
		//create a vector for the map landmark locations predicted within the sensor range of the particle
		vector<LandmarkObs> predicted_landmarks;

		// for each map landmark
		for (int j = 0; j < map_landmarks.landmark_list.size(); j++){

			// get map landmark id, x, y coordinates  (global coordinates)
			double lm_x = map_landmarks.landmark_list[j].x_f;
			double lm_y = map_landmarks.landmark_list[j].y_f;
			int lm_id = map_landmarks.landmark_list[j].id_i;

			// find landmarks within the sensor range of the particle considering a circular region
			//if (dist(lm_x, lm_y, p.x, p.y) <= sensor_range) 

			// find landmarks within the sensor range of the particle considering a rectangular region
			if (fabs(lm_x - p_x) <= sensor_range && fabs(lm_y - p_y) <= sensor_range)
			{
				// create landmark in vehicle coordinate
				LandmarkObs landmark;

				landmark.id = lm_id;
				landmark.x = lm_x;
				landmark.y = lm_y;

				// add predicted to vector predicted_landmarks
				predicted_landmarks.push_back(landmark);
			}
		}
		// create a copy of a list of observations from vehicle coordinates to map coodinates (global coordinates)
		vector<LandmarkObs> transformed_obs;

		for (int k = 0; k < observations.size(); k++){

			//  rotation AND translation transformation transform_x = T * x 
			// using homogeneous transformation matrix T = (cos(theta), -sin(theta), x_t
			//												sin(theta), cos(theta) , y_t
			//														 0, 		  0, 1 )
			//

			// current observation
			LandmarkObs obs = observations[k];

			LandmarkObs obs_landmark;

			// roate and translate
			obs_landmark.x = cos(p_theta) * obs.x - sin(p_theta) * obs.y + p_x;
			obs_landmark.y = sin(p_theta) * obs.x + cos(p_theta) * obs.y + p_y;
			obs_landmark.id = obs.id;

			// add transform landmark to a vector
			transformed_obs.push_back(obs_landmark);
		}

		// perform dataAssociation for predicted_landmarks and transformed_obs on current particle
		dataAssociation(predicted_landmarks, transformed_obs);

		// initialize the weight
		p.weight = 1.0;

		for ( int j = 0; j < transformed_obs.size(); j++){
			// placeholder for the observation and associated prediction coordinates
			double obs_x, obs_y, pred_x, pred_y;

			obs_x = transformed_obs[j].x;
			obs_y = transformed_obs[j].y;

			int associated_pred_id = transformed_obs[j].id;

			// get the x,y coordiantes of the prediction asscciated with the current observation
			for (int k = 0; k < predicted_landmarks.size(); k++){
				if (predicted_landmarks[k].id == associated_pred_id){
					pred_x = predicted_landmarks[k].x; // mean of x position
					pred_y = predicted_landmarks[k].y; // mean of y position
				}
			}

			// calculate the weight for this observation with multi-variant Gaussian 
			double sigma_x = std_landmark[0];
			double sigma_y = std_landmark[1];
			double obs_weight = (1/ (2 * M_PI * sigma_x * sigma_y) ) * exp(- ( pow(pred_x-obs_x,2) / (2 * pow(sigma_x,2)) + pow(pred_y-obs_y,2) / (2 * pow(sigma_y,2)) ) );

			//calculate the total observation weight
			p.weight *= obs_weight;				
		}
		weights[i] = p.weight;
		*/


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	// initialize random engine
	default_random_engine gen;

	vector<Particle> new_particles;

	//sample num_particles from weighted discrete distribution
	discrete_distribution<int> discr_dist(weights.begin(), weights.end());

	// resampling wheel
	for (int i = 0; i < num_particles; i++)
	{
		new_particles.push_back(particles[discr_dist(gen)]);
	}
	particles = new_particles;
}


Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
