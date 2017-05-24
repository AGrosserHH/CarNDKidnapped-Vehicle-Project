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

#include "particle_filter.h"

#define NUMBER_OF_PARTICLES 101

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
       num_particles=NUMBER_OF_PARTICLES;
    
    default_random_engine gen;   
    gen.seed(234);
    
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);   
                
     for (int i=0; i<num_particles;i++){
         Particle p;
     
        p.id=i;
        p.x=dist_x(gen);
        p.y=dist_y(gen);
        p.theta= dist_theta(gen);        
        p.weight=1.0;
        
        particles.push_back(p);   
       
    }
    
    is_initialized=true;
    
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    
    
    default_random_engine gen;   
    gen.seed(234);
    
    const double vel_d_t = velocity * delta_t;
    const double yaw_d_t = yaw_rate * delta_t;
    const double vel_raw = velocity/yaw_rate;
    
    num_particles=NUMBER_OF_PARTICLES;
        
    normal_distribution<double> dist_x(0.0, std_pos[0]);
    normal_distribution<double> dist_y(0.0, std_pos[1]);
    normal_distribution<double> dist_theta(0.0, std_pos[2]);
    
    for (int i=0; i<num_particles;i++){
        
        if (fabs(yaw_rate) < 0.001){
            particles[i].x += vel_d_t * cos(particles[i].theta);
            particles[i].y += vel_d_t * sin(particles[i].theta);
        } else{
        
            particles[i].x += vel_raw * ((sin((particles[i].theta+yaw_d_t)))-(sin(particles[i].theta)));
            particles[i].y += vel_raw * ((cos(particles[i].theta))- (cos((particles[i].theta+yaw_d_t))));
            particles[i].theta += yaw_d_t;
        }
        
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
  
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
        
    for (int i = 0; i<observations.size(); i++){
        
        LandmarkObs o = observations[i];
        
        double min_dist = numeric_limits<double>::max();
        
        int pred_id=-1;
        
        for (int j=0; j<predicted.size(); j++){
            LandmarkObs p = predicted[j];
            
            double cur_dist = sqrt((p.x - o.x) * (p.x - o.x) + (p.y - o.y) * (p.y - o.y));
            
            if (cur_dist<min_dist) {
                min_dist = cur_dist;
                pred_id=p.id;     
                
            }
        }
        
        observations[i].id=pred_id;
    }    
       
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
    
    
    num_particles=NUMBER_OF_PARTICLES;
    for (int i = 0; i < num_particles; i++) {
        
        //get particles coordinates
        double px = particles[i].x;
        double py = particles[i].y;
        double ptheta = particles[i].theta;
                
        // create vector Landmarks predicted
        vector<LandmarkObs> predictions;    
        
        for (int j=0; j<map_landmarks.landmark_list.size(); j++){
            
            int lm_id = map_landmarks.landmark_list[j].id_i ;
	    float lm_x= map_landmarks.landmark_list[j].x_f;
	    float lm_y= map_landmarks.landmark_list[j].y_f;          
            
            if ((fabs(lm_x-px)<=sensor_range) && (fabs(lm_x-px)<=sensor_range)){
                predictions.push_back(LandmarkObs{lm_id, lm_x, lm_y});                
            }            
        }            
        
        vector<LandmarkObs> transformed_obs;
        for (int j=0;j<observations.size(); j++){
            
            int trans_id = observations[j].id; 
            double trans_x = cos(ptheta)*observations[j].x-sin(ptheta)*observations[j].y+px;
            double trans_y = sin(ptheta)*observations[j].x+cos(ptheta)*observations[j].y+py;

            transformed_obs.push_back(LandmarkObs{trans_id, trans_x, trans_y});
        }
        
        dataAssociation(predictions, transformed_obs);
        
        // reinit weight
        particles[i].weight = 1.0;
        
        double o_x, o_y, pr_x, pr_y;
        int asso_pred;
        
        for (int k=0; k < transformed_obs.size(); k++) {
            
            asso_pred = transformed_obs[k].id;
            o_x = transformed_obs[k].x;
            o_y = transformed_obs[k].y;
            
            for (int l=0; l < predictions.size();l++){
                
                if (predictions[l].id==asso_pred) {                    
                    pr_x = predictions[l].x;
                    pr_y = predictions[l].y;                
                    
                }                
                
            }
            
            double std_0 = std_landmark[0];
            double std_1 = std_landmark[1];
            double obs_w = ( 1/(2*M_PI*std_0*std_1)) * exp( -( pow(pr_x-o_x,2)/(2*pow(std_0, 2)) + (pow(pr_y-o_y,2)/(2*pow(std_1, 2))) ) );
            
            particles[i].weight*=obs_w;
            
        }        
    }
    
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
        num_particles=NUMBER_OF_PARTICLES;
        default_random_engine gen;
        vector<Particle> new_particles;
        
        // get weights
        vector<double> weights;
        for (int i = 0; i < num_particles; i++) {            
            weights.push_back(particles[i].weight);            
        }
        
        // generate random distribution for resampling
        uniform_int_distribution<int> unii_dist(0, num_particles-1);
        auto index = unii_dist(gen);

        // get max weight
        double m_weight = *max_element(weights.begin(), weights.end());

        // random distribution (0.0, max_weight)
        uniform_real_distribution<double> unir_dist(0.0, m_weight);
        
        double beta = 0.0;
        
        // spin the resample wheel!
        for (int i = 0; i < num_particles; i++) {
          beta += unir_dist(gen) * 2.0;
          while (beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
          }
          new_particles.push_back(particles[index]);
        }

        particles = new_particles;      
         
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
        num_particles=NUMBER_OF_PARTICLES;
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
