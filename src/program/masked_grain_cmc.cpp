//----------------------------------------------------------------------------------------------------
// Program for calculating individual grain properties in HAMR simulations via masked Monte Carlo
//----------------------------------------------------------------------------------------------------

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "grains.hpp"

namespace program{

// Main function to distinguish grains via use of a mask
void masked_grain_cmc(){

   // Vector ensuring all grains in system are constrained
	std::vector<bool> constrained(grains::num_grains, true);

   // Size of Phi and Theta grain information vector
	double phi_theta_config=2*grains::num_grains;

   // Initialization of Phi and Theta vector for each grain
	std::vector<double> phi_theta_constraints(phi_theta_config, 0.0);

   // Allocate constraint angles to each grain according to user input values
	for (int grain_parameter=0; grain_parameter<=phi_theta_config; grain_parameter=grain_parameter+2){

		phi_theta_constraints[grain_parameter] = sim::constraint_phi;
		phi_theta_constraints[grain_parameter+1] = sim::constraint_theta;

	}

   // Execute masked Monte Carlo integrator based on grain mask
	montecarlo::initialise_masked_cmc_mc(grains::num_grains, atoms::grain_array, constrained, phi_theta_constraints); 

	return;
    
}

} // End of 'program' namespace