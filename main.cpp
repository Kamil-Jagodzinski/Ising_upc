#include "utils.h"

void runProgram(int rank, int num_procs, int grid_size, double J,
                double B, long long iterations, long long repeat);

int main(int argc, char* argv[]) {
    upcxx::init();

    int rank = upcxx::rank_me();
    int num_procs = upcxx::rank_n();

    int grid_size;
    double J, B;
    long long iterations, repeat;
    bool stay_in_GUI = true;
    readParametersFromFile(grid_size, J, B, iterations, repeat);

    if (rank == 0) {
        while (true && stay_in_GUI) {
            std::cout << "====================================================" << std::endl;
            std::cout << "||         SELECT ACTION FROM LIST BELOW:         ||" << std::endl;
            std::cout << "|| 1. Change grid_size      4. Change iterations  ||" << std::endl;
            std::cout << "|| 2.     Change J          5.   Change repeat    ||" << std::endl;
            std::cout << "|| 3.     Change B          6.    Run program     ||" << std::endl;
            std::cout << "||                                                ||" << std::endl;
            std::cout << "||                                                ||" << std::endl;
            std::cout << "====================================================" << std::endl;
            std::cout << "                Currently set parameters:" << std::endl;
            std::cout << "grid_size = " << grid_size << ",\t J = " << J << ",\t B = " << B
                      << ",\niterations = " << iterations << ",\nrepeat = " << repeat << std::endl;

            std::cout << "\n>>  Choose option: \n>>";

            int option;
            std::cin >> option;

            switch (option) {
                case 1:
                    system("clear");
                    std::cout << "Enter new Net Size: ";
                    while (!(std::cin >> grid_size) || grid_size <= 0) {
                        std::cout << "Invalid input. Please enter a positive integer for Net Size: ";
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    saveParametersToFile(grid_size, J, B, iterations, repeat);
                    break;
                case 2:
                    system("clear");
                    std::cout << "Enter new J (-1 to 1): ";
                    while (!(std::cin >> J) || J < -1 || J > 1) {
                        std::cout << "Invalid input. Please enter a value between -1 and 1 for J: ";
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    saveParametersToFile(grid_size, J, B, iterations, repeat);
                    break;
                case 3:
                    system("clear");
                    std::cout << "Enter new B (-1 to 1): ";
                    while (!(std::cin >> B) || B < -1 || B > 1) {
                        std::cout << "Invalid input. Please enter a value between -1 and 1 for B: ";
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    saveParametersToFile(grid_size, J, B, iterations, repeat);
                    break;
                case 4:
                    system("clear");
                    std::cout << "Enter new iterations: ";
                    while (!(std::cin >> iterations) || iterations <= 0) {
                        std::cout << "Invalid input. Please enter a positive integer for iterations: ";
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    saveParametersToFile(grid_size, J, B, iterations, repeat);
                    break;
                case 5:
                    system("clear");
                    std::cout << "Enter new repeat: ";
                    while (!(std::cin >> repeat) || repeat <= 0) {
                        std::cout << "Invalid input. Please enter a positive integer for repeat: ";
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    saveParametersToFile(grid_size, J, B, iterations, repeat);
                    break;
                case 6:
                    system("clear");
                    std::cout << "Running the program" << std::endl;
                    stay_in_GUI = false;
                    break;
                default:
                    std::cout << "Invalid option" << std::endl;
                    break;
            }
        }
    }

    upcxx::barrier();

    runProgram(rank, num_procs, grid_size, J, B, iterations, repeat);

    upcxx::finalize();
    return 0;
}

void runProgram(int rank, int num_procs, int grid_size, double J,
                double B, long long iterations, long long repeat) {
    int row_size = grid_size;
    int iters = iterations;
    int rows_per_proc = row_size / num_procs;
    std::string dir_name;
    std::mt19937 gen(rank);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::uniform_int_distribution<int> disInt(0, INT_MAX);

    for (int rep = 0; rep < repeat; rep++) {
        if (rank == 0) {
            dir_name = createFolderWithTimestampName(rep);
            if (dir_name == "ERROR") {
                std::cout << "ERROR during creating dir" << std::endl;
                return;
            }
        }

        upcxx::global_ptr<int> cluster = upcxx::new_array<int>(row_size * rows_per_proc);
        upcxx::global_ptr<int> recv_buffer = upcxx::new_array<int>(row_size * row_size);

        upcxx::barrier();

        upcxx::broadcast(&J, 1, 0).wait();
        upcxx::broadcast(&B, 1, 0).wait();
        upcxx::broadcast(&grid_size, 1, 0).wait();
        upcxx::broadcast(&iterations, 1, 0).wait();
        upcxx::broadcast(&repeat, 1, 0).wait();

        for (int i = 0; i < iters + num_procs; i += num_procs) {
            int idx = disInt(gen) % (rows_per_proc * row_size) + rank * rows_per_proc * row_size;

            double delta = upcxx::rpc(rank, calculateEnergyChange, recv_buffer, idx, row_size,
                                      rows_per_proc, num_procs).wait();
            double p = 0.0;

            if (delta < 0.0) {
                p = 1.0;
            } else {
                p = std::exp(-delta / J);
            }

            double rnd = dis(gen);

            if (rnd < p) {
                upcxx::rpc(rank, flipSpin, recv_buffer, idx, row_size,
                           rows_per_proc, num_procs).wait();
            }

            if (i % num_procs == 0) {
                upcxx::barrier();
                if (rank == 0) {
                    saveSpinToFile(recv_buffer, row_size, i / num_procs, dir_name);
                }
            }
        }

        upcxx::delete_array(cluster);
        upcxx::delete_array(recv_buffer);

        upcxx::barrier();
    }
}