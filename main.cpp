#include "utils.h"
#include <upcxx/upcxx.hpp>
#include <random>
#include <cmath>
#include <climits>
#include <algorithm>

int sys_retval;

void runProgram(int rank, int num_procs, int grid_size, double J, 
                double B, long long iterations, long long repeat);

int main(int argc, char *argv[]) {
    upcxx::init();

    int rank = upcxx::rank_me();
    int num_procs = upcxx::rank_n();
    
    int grid_size; 
    double J, B; 
    long long iterations, repeat;
    bool stay_in_GUI = true;

    // Assuming this function is adjusted to work with UPC++ as necessary
    readParametersFromFile( grid_size, J, B, iterations, repeat);

    // Ensure all other processes wait until GUI session is done
    upcxx::barrier();

    // Run the computations
    runProgram(rank, num_procs, grid_size, J, B, iterations, repeat);
    
    // Finalize program
    upcxx::finalize();
    return 0;
}



void runProgram(int rank, int num_procs, int grid_size, double J, double B, long long iterations, long long repeat) {

    
    // Initialization
    int row_size = grid_size;
    int iters = iterations;
    int rows_per_proc = row_size/num_procs;
    std::string dir_name;
    std::mt19937 gen(rank); 
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::uniform_int_distribution<int> disInt(0, INT_MAX);
    upcxx::global_ptr<int> myGlobalPtr = nullptr;

    

    // Using UPC++ dist_object to communicate the variables between ranks
    upcxx::dist_object<int> d_grid_size(grid_size);
    upcxx::dist_object<double> d_J(J);
    upcxx::dist_object<double> d_B(B);
    upcxx::dist_object<long long> d_iterations(iterations);
    upcxx::dist_object<long long> d_repeat(repeat);

    
    for(int rep=0; rep<repeat; rep++){
        
        if(rank == 0){
            dir_name = createFolderWithTimestampName(rep);
            if ( dir_name == "ERROR" ){
               std::cout << "ERROR during creating dir" << std::endl;
                return ;
            }
        }
        
        if(rank == 0){
            myGlobalPtr = upcxx::new_array<int>(grid_size * grid_size);
            for(int i=0; i<grid_size; i++){
                for(int j=0; j<grid_size; j++){
                    upcxx::rput(1, myGlobalPtr + (i * row_size) + j).wait();
                }
            }
        }

        // Broadcast using upcxx::broadcast and wait for it to complete
        myGlobalPtr = upcxx::broadcast(myGlobalPtr, 0).wait();
        J = upcxx::broadcast(*d_J, 0).wait();
        B = upcxx::broadcast(*d_B, 0).wait();
        grid_size = upcxx::broadcast(*d_grid_size, 0).wait();
        iterations = upcxx::broadcast(*d_iterations, 0).wait();
        repeat = upcxx::broadcast(*d_repeat, 0).wait();

        upcxx::barrier();


        // Run the simulation
        for(int i=0; i<iters+num_procs; i+=num_procs) {
            int idx = disInt(gen) % (rows_per_proc * row_size) + rank * rows_per_proc * row_size;

            double delta = calculateEnergyChange(myGlobalPtr, idx, row_size, rows_per_proc, num_procs).wait();
            double p = (delta < 0.0) ? 1.0 : exp(-delta);


            if( dis(gen) < p ){
                (upcxx::rget(myGlobalPtr+idx).wait() == 0) ? upcxx::rput(1, myGlobalPtr + idx).wait() : upcxx::rput(0, myGlobalPtr + idx).wait();
            }
            upcxx::barrier();
            

            // Save data
            if(rank == 0 && (i == iters-1 || i%(iters/10) < num_procs)) {
                saveGrid(myGlobalPtr, row_size, dir_name);
            }
            upcxx::barrier();


            if(rank == 0 && (i == iters-1 || i%(iters/100) < num_procs)) {
                saveEnergy(energy(myGlobalPtr, J, B, row_size), dir_name);
                saveMag(avgMagnetism(myGlobalPtr, row_size * rows_per_proc * num_procs), dir_name);
            }

            upcxx::barrier();
        }

        
    }

    // Cleanup
    upcxx::barrier();
    if(rank == 0){
        upcxx::delete_(myGlobalPtr);
        std::cout << "Symulacja zakończona" << std::endl;
    }
}
