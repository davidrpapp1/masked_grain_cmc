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

void masked_grain_cmc(){

	std::vector<bool> constrained(grains::num_grains, true);

	std::cout << grains::num_grains << std::endl;

	double phi_theta_config=2*grains::num_grains;

	std::vector<double> phi_theta_constraints(phi_theta_config, 0.0);

	for (int grain_parameter=0; grain_parameter<=phi_theta_config; grain_parameter=grain_parameter+2){

		phi_theta_constraints[grain_parameter] = sim::constraint_phi;
		phi_theta_constraints[grain_parameter+1] = sim::constraint_theta;

	}

	montecarlo::initialise_masked_cmc_mc(grains::num_grains, atoms::grain_array, constrained, phi_theta_constraints); 

	return;
    
}

}