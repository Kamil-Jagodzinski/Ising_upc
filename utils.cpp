#include "utils.h"

upcxx::global_ptr<int> initializeGrid(int grid_size) {
    upcxx::global_ptr<int> grid = nullptr;
    if (upcxx::rank_me() == 0) {
        // Inicjalizacja siatki spinów tylko na procesie o ranku 0
        grid = upcxx::new_array<int>(grid_size * grid_size);
        int* local_grid = grid.local();
        
        // Wypełnienie siatki spinów jedynkami
        for (int i = 0; i < grid_size * grid_size; i++) {
            local_grid[i] = 1;
        }
    }
    
    // Broadcast siatki spinów do wszystkich procesów
    grid = upcxx::broadcast(grid, 0).wait();
    
    return grid;
}

void printVector2D(const upcxx::global_ptr<int>& grid, int grid_size) {
    if (upcxx::rank_me() == 0) {
        int* local_grid = grid.local();
        std::cout << "Model grid:" << std::endl;
        for (int i = 0; i < grid_size; i++) {
            for (int j = 0; j < grid_size; j++) {
                std::cout << local_grid[i * grid_size + j] << " ";
            }
            std::cout << std::endl;
        }
    }
    upcxx::barrier();
}

double calculateEnergyChange(const upcxx::global_ptr<int>& grid, int idx, int row_size, int rows_per_proc, int num_procs) {
    int left_idx = (idx % row_size != 0) ? idx - 1 : idx + row_size - 1;
    int right_idx = (idx % row_size != row_size - 1) ? idx + 1 : idx - row_size + 1;
    int up_idx = (idx >= row_size) ? idx - row_size : idx + rows_per_proc * row_size - row_size;
    int down_idx = (idx < rows_per_proc * row_size - row_size) ? idx + row_size : idx % row_size;

    int spin = upcxx::rget(grid + idx).wait() == 1 ? 1 : -1;
    int left_spin = upcxx::rget(grid + left_idx).wait() == 1 ? 1 : -1;
    int right_spin = upcxx::rget(grid + right_idx).wait() == 1 ? 1 : -1;
    int up_spin = upcxx::rget(grid + up_idx).wait() == 1 ? 1 : -1;
    int down_spin = upcxx::rget(grid + down_idx).wait() == 1 ? 1 : -1;

    return 2.0 * spin * (left_spin + right_spin + up_spin + down_spin);
}

double single_spin_energy(int index, const upcxx::global_ptr<int>& grid, int row_size, double J, double B) {
    int energyNeigh = -2; // Initial value of -2 for binary grid values (0 and 1)

    int x = index / row_size;
    int y = index % row_size;

    // Periodic boundary conditions
    int left = (y == 0) ? x * row_size + row_size - 1 : x * row_size + y - 1;
    int right = (y == row_size - 1) ? x * row_size : x * row_size + y + 1;
    int up = (x == 0) ? (row_size - 1) * row_size + y : (x - 1) * row_size + y;
    int down = (x == row_size - 1) ? y : (x + 1) * row_size + y;

    upcxx::future<int> left_future = upcxx::rget(grid + left);
    upcxx::future<int> right_future = upcxx::rget(grid + right);
    upcxx::future<int> up_future = upcxx::rget(grid + up);
    upcxx::future<int> down_future = upcxx::rget(grid + down);

    energyNeigh += left_future.wait() + right_future.wait() + up_future.wait() + down_future.wait();

    return (J * static_cast<double>(upcxx::rget(grid + index).wait()) * static_cast<double>(energyNeigh)
            + B * 0.25 * (upcxx::rget(grid + index).wait() ? 1.0 : -1.0));
}

double energy(const upcxx::global_ptr<int>& grid, double J, double B, int row_size) {
    double sum = 0.0;

    for (int i = 0; i < row_size * row_size; i++) {
        sum += single_spin_energy(i, grid, row_size, J, B);
    }

    return sum;
}

void flipSpin(upcxx::global_ptr<int>& grid, int idx) {
    int* local_grid = grid.local();
    local_grid[idx] = (local_grid[idx] == 0) ? 1 : 0;
}

double avgMagnetism(const upcxx::global_ptr<int>& spinArray, int spinArraySize) {
    double sum = 0.0;

    upcxx::rput(0.0, sum).wait();
    upcxx::barrier();

    for (int i = 0; i < spinArraySize; i++) {
        int spin = upcxx::rget(spinArray + i).wait();
        sum += spin;
    }

    double global_sum = upcxx::reduce_one(sum, upcxx::op_fast_add, 0).wait();
    double avg = global_sum / spinArraySize;

    return avg;
}

void saveGrid(const upcxx::global_ptr<int>& grid, int row_size, std::string folderName) {
    char filename[256];
    const char* cstr = folderName.c_str();
    sprintf(filename, "%s/spins.txt", cstr);
    upcxx::global_ptr<FILE> fp = nullptr;

    if (upcxx::rank_me() == 0) {
        fp = fopen(filename, "a");
        if (fp == nullptr) {
            printf("Error: could not open file for writing.\n");
            return;
        }
    }

    upcxx::barrier();

    for (int i = 0; i < row_size; i++) {
        for (int j = 0; j < row_size; j++) {
            int spin = upcxx::rget(grid + i * row_size + j).wait();

            if (upcxx::rank_me() == 0) {
                fprintf(fp.local(), "%d ", spin);
            }
        }

        if (upcxx::rank_me() == 0) {
            fprintf(fp.local(), "\n");
        }
    }

    upcxx::barrier();

    if (upcxx::rank_me() == 0) {
        fclose(fp).wait();
    }
}

void saveMag(double mg, std::string folderName) {
    char filename[256];
    const char *cstr = folderName.c_str();
    sprintf(filename, "%s/avgMagnetism.txt", cstr);
    FILE* fp = fopen(filename, "a");
    if (fp == NULL) {
        printf("Error: could not open file for writing.\n");
        return;
    }
    fprintf(fp, "%f\n", mg);
    fclose(fp);
}

void saveEnergy(double energy, std::string folderName) {
    char filename[256];
    const char *cstr = folderName.c_str();
    sprintf(filename, "%s/energy.txt", cstr);
    FILE* fp = fopen(filename, "a");
    if (fp == NULL) {
        printf("Error: could not open file for writing.\n");
        return;
    }
    fprintf(fp, "%f\n", energy);
    fclose(fp);
}

std::string createFolderWithTimestampName(int rep){
    // Uzyskaj aktualny czas
    auto currentTime = std::chrono::system_clock::now();
    std::time_t currentTime_t = std::chrono::system_clock::to_time_t(currentTime);

    // Sformatuj czas jako string
    char timestamp[128];
    std::strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", std::localtime(&currentTime_t));
    std::string timestampStr(timestamp);

    // Stwórz ścieżkę do nowego folderu
    std::string folderName =  "result/" + timestampStr + "_" + std::to_string(rep);
    std::filesystem::path folderPath(folderName);

    // Sprawdź czy folder już istnieje
    if (std::filesystem::exists(folderPath))
    {
        std::cerr << "Folder o nazwie " << folderName << " już istnieje.\n";
        return folderName;
    }

    // Stwórz nowy folder
    if (!std::filesystem::create_directory(folderPath))
    {
        std::cerr << "Nie udało się utworzyć folderu " << folderName << "\n";
        return "ERROR";
    }

    std::cout << "Utworzono folder " << folderName << "\n";
    return folderName;
}

void saveParametersToFile(int netSize, double J, double B, long long iters, long long repeat) {
    std::ofstream file("parameters.txt");
    if (!file) {
        std::cout << "Failed to open file for writing." << std::endl;
        return;
    }

    file << "Net Size: " << netSize << std::endl;
    file << "J: " << J << std::endl;
    file << "B: " << B << std::endl;
    file << "Number of iterations: " << iters << std::endl;
    file << "Number repeats: " << repeat << std::endl;

    file.close();
    std::cout << "Parameters saved to file successfully." << std::endl;
}

void readParametersFromFile(int& netSize, double& J, double& B, long long& iters, long long& repeat) {
    std::ifstream file("parameters.txt");
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (line.find("Net Size:") != std::string::npos) {
                netSize = std::stoi(line.substr(line.find(":")+1));
            } else if (line.find("J:") != std::string::npos) {
                J = std::stod(line.substr(line.find(":")+1));
            } else if (line.find("B:") != std::string::npos) {
                B = std::stod(line.substr(line.find(":")+1));
            } else if (line.find("Number of iterations:") != std::string::npos) {
                iters = std::stoll(line.substr(line.find(":")+1));
            } else if (line.find("Number repeats:") != std::string::npos) {
                repeat = std::stoll(line.substr(line.find(":")+1));
            }
        }
        file.close();
    } else {
        std::cerr << "Error: could not open file for reading." << std::endl;
    }
}