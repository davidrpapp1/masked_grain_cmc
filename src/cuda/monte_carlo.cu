
#include "atoms.hpp"
#include "exchange.hpp"
#include "vio.hpp"
#include "cuda_utils.hpp"
#include "internal.hpp"
#include "data.hpp"
#include "material.hpp"

#include "spin_fields.hpp"

#include "monte_carlo.hpp"

#include <vector>

namespace vcuda
{


   //--------------------------------------------------------------------------
   // Function to perform a single Monte Carlo step
   //--------------------------------------------------------------------------
   void mc_step(){

      #ifdef CUDA
        // check for cuda initialization, and initialize if necessary
        if (!internal::mc::initialised) internal::mc::initialise();
        // perform a single monte carlo step
        internal::mc::__mc_step();

      #endif

      return;
   }


    namespace internal
    {
        namespace mc
        {

            bool initialised(false);

            // Number of sublattices
            int M;

            cu_real_t step_size;

            // Array to store the colour of each atom
            std::vector<int> atom_colour;

            // Create a set of lists giving which atoms have which colours
            std::vector< std::vector<int> > colour_list;

            //Index starting point in continuous array
            std::vector<int> sl_start;

            // device array for atoms in (sub lattice)
            int * d_sl_atoms;

            int * d_accepted;

            int64_t             seed;
            curandGenerator_t   gen;

            // Arrays for random spin vector and acceptance probability
            cu_real_t *d_rand_spin;
            cu_real_t *d_rand_accept;

            int colour_split()
            {

                std::cout << "Starting graph colouring" << std::endl;

                int N = ::atoms::num_atoms;

                // Array to store the colour of each atom
                // All are initially uncoloured (=-1)
                std::vector<int> atom_colour;
                atom_colour.assign(N, -1);

                // initial number of colours
                M = 2;

                int Ncoloured = 0;

                // Starting atom
                int iatom = 0;

                // Array to store the number of each colour that nbr the atoms
                std::vector<int> nbr_colours(M);

                do{

                    //std::cout << Ncoloured << "  " << iatom << "  " << atom_colour[iatom] << "  ";
                    // temporay constants for loop start and end indices
                    const int start = ::atoms::neighbour_list_start_index[iatom];
                    const int end   = ::atoms::neighbour_list_end_index[iatom]+1;


                    nbr_colours.assign(M, 0);
                    // loop over all neighbours and count how many of each colour
                    for(int nn = start; nn < end; ++nn){
                        const int atomj = ::atoms::neighbour_list_array[nn]; // get neighbouring atom number
                        int ac = atom_colour[atomj];
                        if (ac >= 0)
                            nbr_colours[ac]++;
                    }

                    // Go backwards through the number of nbring colours
                    // This atoms colour will be the lowest zero colour
                    int lowest = M;
                    for (int i = M-1; i >= 0; i--)
                        if (nbr_colours[i] == 0)
                            lowest = i;


                    // If all there is no available colour increase the total and give the atom the new colour
                    if (lowest == M)
                        M++;

                    atom_colour[iatom] = lowest;
                    Ncoloured++;
                    //std::cout << "lowest = " << lowest << std::endl;


                    // Select the next atom from the next uncoloured neighbour
                    int next_atom = -1;
                    // loop over all neighbours
                    for(int nn = start; nn < end; ++nn){
                        const int atomj = ::atoms::neighbour_list_array[nn]; // get neighbouring atom number
                        int ac = atom_colour[atomj];
                        // find the lowest uncoloured nbr
                        if (ac < 0){
                            next_atom = atomj;
                            break;
                        }
                    }

                    //if all neighbours are coloured select from list
                    if (next_atom == -1){
                        for ( int i = 0; i < N; i++){
                            int ac = atom_colour[i];
                            // find the lowest uncoloured nbr
                            if (ac < 0){
                                next_atom = i;
                                break;
                            }
                        }
                    }


                    // Move onto the next selected atom
                    iatom = next_atom;

                    // If we can't find another uncoloured atom, end the do while
                    if (next_atom == -1){
                        std::cout << "No more uncoloured atoms. Fraction coloured = " << Ncoloured << " / " << N << std::endl;
                        break;
                    }


                }while(Ncoloured < N);


                // Create a set of lists giving which atoms have which colours
                colour_list.resize(M);

                for ( int i = 0; i < N; i++){
                    int ac = atom_colour[i];
                    colour_list[ac].push_back(i);
                }

                for (int i = 0; i < M; i++){
                    std::cout << "Colour: " << i << ", natoms = " << colour_list[i].size() << std::endl;
                }


                return 0;
            }

            int initialise()
            {

                curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
                curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);


                cudaMalloc((void**)&d_rand_spin, 3*::atoms::num_atoms * sizeof(cu_real_t));
                cudaMalloc((void**)&d_rand_accept, 3*::atoms::num_atoms * sizeof(cu_real_t));

                cudaMalloc((void**)&d_sl_atoms, ::atoms::num_atoms * sizeof(int));
                cudaMalloc((void**)&d_accepted, ::atoms::num_atoms * sizeof(int));

                colour_split();

                std::vector<int> h_sl_atoms;
                h_sl_atoms.resize(::atoms::num_atoms);
                int idx = 0;
                for ( int i = 0; i < M; i++){
                    sl_start.push_back(idx);
                    for ( int j = 0; j < colour_list[i].size(); j++){
                        h_sl_atoms[idx] = colour_list[i][j];
                        idx++;
                    }
                }

                cudaMemcpy(d_sl_atoms, h_sl_atoms.data(), ::atoms::num_atoms * sizeof(int), cudaMemcpyHostToDevice);

                std::cout << "Trying a step..."<< std::endl;

                __mc_step();
                std::cout << "Done"<< std::endl;

                initialised = true;

                return 0;
            }

            void finalise()
            {

                curandDestroyGenerator(gen);

                cudaFree(d_rand_spin);
                cudaFree(d_rand_accept);

                cudaFree(d_sl_atoms);
                cudaFree(d_accepted);
            }

            __global__ void monte_carlo_sublattice_step(
                    const int sl_start,
                    int * sl_atoms,
                    int * material,
                    vcuda::internal::material_parameters_t * material_params,
                    cu_real_t *rand_spin,
                    cu_real_t *rand_accept,
                    int * accepted,
                    cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
                    cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
                    const cu_real_t step_size, const cu_real_t global_temperature, const size_t N){

                // Loop over blocks for large systems > ~100k spins
                for ( size_t i = blockIdx.x * blockDim.x + threadIdx.x;
                        i < N;
                        i += blockDim.x * gridDim.x)
                {

                    // load an atom from the current sublattice
                    int atom = sl_atoms[i+sl_start];

                    // Get material of atom
                    int mid = material[atom];

                    // Load parameters to local variables from memory
                    cu::material_parameters_t mat = material_params[mid];

                    // load spin direction to registers for later multiple reuse
                    cu_real_t sx = x_spin[atom];
                    cu_real_t sy = y_spin[atom];
                    cu_real_t sz = z_spin[atom];

                    // new spin direction
                    cu_real_t nsx, nsy, nsz;

                    // run chosen move type
                    // Select algorithm using case statement

                    
                    switch(montecarlo::algorithm){

                        case adaptive:
                        {      
                            
                            nsx = rand_spin[atom];      // old_spin[0] + mtrandom::gaussian() * montecarlo::internal::adaptive_sigma;
                            nsy = rand_spin[atom+N];    // old_spin[1] + mtrandom::gaussian() * montecarlo::internal::adaptive_sigma;
                            nsz = rand_spin[atom+2*N];  // old_spin[2] + mtrandom::gaussian() * montecarlo::internal::adaptive_sigma;
            
                            cu_real_t mod_s = sqrt(nsx*nsx + nsy*nsy + nsz*nsz);
            
                            nsx /= mod_s;
                            nsy /= mod_s;
                            nsz /= mod_s;
                            break;
                        }
                        case spin_flip:
                            nsx = -sx;
                            nsy = -sy;
                            nsz = -sz;
                            break;
                        
                        case uniform:
                        {      
                            nsx = rand_spin[atom];
                            nsy = rand_spin[atom+N];
                            nsz = rand_spin[atom+2*N];
            
                            cu_real_t mod_s = sqrt(nsx*nsx + nsy*nsy + nsz*nsz);
            
                            nsx /= mod_s;
                            nsy /= mod_s;
                            nsz /= mod_s;
                            break;
                        }
                        
                        case angle:
                        {      
                            
                            nsx = rand_spin[atom];      // old_spin[0] + mtrandom::gaussian() * montecarlo::internal::delta_angle;
                            nsy = rand_spin[atom+N];    // old_spin[1] + mtrandom::gaussian() * montecarlo::internal::delta_angle;
                            nsz = rand_spin[atom+2*N];  // old_spin[2] + mtrandom::gaussian() * montecarlo::internal::delta_angle;
            
                            cu_real_t mod_s = sqrt(nsx*nsx + nsy*nsy + nsz*nsz);
            
                            nsx /= mod_s;
                            nsy /= mod_s;
                            nsz /= mod_s;
                            break;
                        }
                        case hinzke_nowak:
                        {
                            switch(pick_move){
                                case 0: // spin flip
                                    nsx = -sx;
                                    nsy = -sy;
                                    nsz = -sz;
                                    break;
                                case 1: // uniform
                                {      
                                    nsx = rand_spin[atom];
                                    nsy = rand_spin[atom+N];
                                    nsz = rand_spin[atom+2*N];
                    
                                    cu_real_t mod_s = sqrt(nsx*nsx + nsy*nsy + nsz*nsz);
                    
                                    nsx /= mod_s;
                                    nsy /= mod_s;
                                    nsz /= mod_s;
                                    break;
                                }
                                
                                case 2: // angle 
                                {      
                            
                                    nsx = rand_spin[atom];      // old_spin[0] + mtrandom::gaussian() * montecarlo::internal::delta_angle;
                                    nsy = rand_spin[atom+N];    // old_spin[1] + mtrandom::gaussian() * montecarlo::internal::delta_angle;
                                    nsz = rand_spin[atom+2*N];  // old_spin[2] + mtrandom::gaussian() * montecarlo::internal::delta_angle;
                    
                                    cu_real_t mod_s = sqrt(nsx*nsx + nsy*nsy + nsz*nsz);
                    
                                    nsx /= mod_s;
                                    nsy /= mod_s;
                                    nsz /= mod_s;
                                    break;
                                }
                                default: // angle
                                {      
                            
                                    nsx = rand_spin[atom];      // old_spin[0] + mtrandom::gaussian() * montecarlo::internal::delta_angle;
                                    nsy = rand_spin[atom+N];    // old_spin[1] + mtrandom::gaussian() * montecarlo::internal::delta_angle;
                                    nsz = rand_spin[atom+2*N];  // old_spin[2] + mtrandom::gaussian() * montecarlo::internal::delta_angle;
                    
                                    cu_real_t mod_s = sqrt(nsx*nsx + nsy*nsy + nsz*nsz);
                    
                                    nsx /= mod_s;
                                    nsy /= mod_s;
                                    nsz /= mod_s;
                                    break;
                                }
                            }
                            break;
                        }
                        default: // adaptive
                        {      
                            
                            nsx = rand_spin[atom];      // old_spin[0] + mtrandom::gaussian() * montecarlo::internal::adaptive_sigma;
                            nsy = rand_spin[atom+N];    // old_spin[1] + mtrandom::gaussian() * montecarlo::internal::adaptive_sigma;
                            nsz = rand_spin[atom+2*N];  // old_spin[2] + mtrandom::gaussian() * montecarlo::internal::adaptive_sigma;
            
                            cu_real_t mod_s = sqrt(nsx*nsx + nsy*nsy + nsz*nsz);
            
                            nsx /= mod_s;
                            nsy /= mod_s;
                            nsz /= mod_s;
                            break;
                        }
                    }

                    // Calculate current energy
                    cu_real_t Eold = ::vcuda::internal::uniaxial_anisotropy_energy(mat, sx, sy, sz);
                    cu_real_t Enew = ::vcuda::internal::uniaxial_anisotropy_energy(mat, nsx, nsy, nsz);

                    cu_real_t dE = (Enew - Eold)/(1.38064852e-23*global_temperature);

                    cu_real_t r_accept = rand_accept[atom];
                    if ( r_accept < exp(-dE) ){
                        x_spin[atom] = nsx;
                        y_spin[atom] = nsy;
                        z_spin[atom] = nsz;
                        accepted[i + sl_start] = 1;
                    }
                }
            }

            void __mc_step()
            {
                // Check for cuda errors in file, line
                check_cuda_errors (__FILE__, __LINE__);

                // generate 3 random doubles per atom
                curandGenerateNormalDouble( gen, d_rand_spin, 3*::atoms::num_atoms, 0.0, 1.0);
                curandGenerateUniformDouble( gen, d_rand_accept, 3*::atoms::num_atoms);

                cudaMemset(d_accepted, 0, ::atoms::num_atoms*sizeof(int));

                // Calculate external fields (fixed for integration step)
                cu::update_external_fields ();

                //Iterate over all the sublattices
                for( int i = 0; i < M; i++){

                    monte_carlo_sublattice_step <<< cu::grid_size, cu::block_size >>> (
                            sl_start[i],
                            d_sl_atoms,
                            ::cu::atoms::d_materials, cu::mp::d_material_params,
                            d_rand_spin, d_rand_accept,
                            d_accepted,
                            ::cu::atoms::d_x_spin, ::cu::atoms::d_y_spin, ::cu::atoms::d_z_spin,
                            ::cu::d_x_external_field, ::cu::d_y_external_field, ::cu::d_z_external_field,
                            step_size, sim::temperature, colour_list[i].size());




                }
            }
        }
    }
}
