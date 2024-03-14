#include<cmath>
#include<vector>
#include <iostream>


// Universal constants G and speed of light c

const double G = 1;

const double c = 1; 

const double dt = 1; 

const int num_of_smarticles = 200; 

//kpc
const float max_radial_extent = 30;

//number of rings 
const int ring_num = 5;



// trackable parameters of simulation particles 
std::vector<std::vector<double>> particle_positions;

std::vector<std::vector<double>> particle_velocities;

std::vector<double> particle_accelerations;

std::vector<float> particle_masses;


// trackable parameters of the blackhole 

const float BH_mass = 100; 

std::vector<float> black_hole_position(2);

// we're going to have to calculate the schwarzchild radius of our non-rotating Schwarzschild blackhole 
// this radius defines the limit at wich the forces that push outwards are overpowered by gravity 

double calculate_rs(float BH_m,double G, double c){

    return ((2*G*BH_m)/pow(c,2));
}

// first the initial conditions need to be created 
// we'll begin with particles arranged in prgressively larger circular rings 
// the rings themselves will be at constant distances to eachother with the BH at the center

std::vector<std::vector<double>> initialize_positions(int number_of_particles, float r_max, int r_nums){

    std::vector<std::vector<double>> initial_pos_array(number_of_particles, std::vector<double>(2));

    std::vector<double> radial_distances;

    int r_increment = r_max/r_nums;

    for (double r = 0; r < r_max; r += r_increment){

        radial_distances.push_back(r);
    }

    int number_of_rings = radial_distances.size();

    // we would like approximately the same number of particles per ring

    int number_of_particles_per_ring = round(number_of_particles/radial_distances.size());

    int particles_left_over = (number_of_particles_per_ring*number_of_rings) - number_of_particles;

    // 2pi radians = 360 degrees 

    float angle_increment = 2 * M_PI / number_of_particles_per_ring;

    int particle_num = 0; 

    for (int i = 0; i < number_of_rings; ++i){

        float at_angle = 0;

        // the radial distances are arranged in ascending order 
        // so we iterate outwards from the innermost radial distance 

        float at_radius = radial_distances[i];


        int particles_in_this_ring = number_of_particles_per_ring;

        // we're adding leftover particles to the innermost ring 
        if (i == 1)
        {
            particles_in_this_ring += particles_left_over;
        }
        
        for(int p; p < particles_in_this_ring; p++){
            
            // x coordinate
            initial_pos_array[particle_num][0] = at_radius * std::cos(at_angle);
            // y coordinate 
            initial_pos_array[particle_num][1] = at_radius * std::sin(at_angle);

            ++particle_num;

            at_angle += angle_increment;


        } // terminating loop iterating over particles in the ring

    } // terminating loop iterating over rings 

    
    return initial_pos_array;

}   



std::vector<std::vector<double>> initialize_velocities(int number_of_particles){

     std::vector<std::vector<double>> initial_vel(number_of_particles, std::vector<double>(2));

     for (int i = 0; i < number_of_particles; i++){

        double x_comp_v = 0.0;
        double y_comp_v = 0.0; 

        initial_vel[i][0] = x_comp_v;
        initial_vel[i][1] = y_comp_v;

     }

    return initial_vel;

}

int main(){

    // size or radial extent of blackhole 
    const double r_schwarzschild = calculate_rs(BH_mass,G,c);  
    
    // start by placing the blackhole at the center
    black_hole_position[0]=0;
    black_hole_position[1]=0;

    particle_positions = initialize_positions(num_of_smarticles, max_radial_extent, ring_num);

    std::vector<std::vector<double>> initial_vel(num_of_smarticles, std::vector<double>(2));

    particle_velocities = initialize_velocities(num_of_smarticles);

    particle_accelerations;

    return 0;

    

}// terminating main function 







