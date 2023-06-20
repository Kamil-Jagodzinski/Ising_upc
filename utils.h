#pragma once
#include <upcxx/upcxx.hpp>
#include <iostream>
#include <random>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib> 
#include <bits/stdc++.h>


upcxx::global_ptr<int> initializeGrid(int grid_size);

void printVector2D(const upcxx::global_ptr<int>& grid, int grid_size);

double calculateEnergyChange(const upcxx::global_ptr<int>& grid, int idx, int row_size, int rows_per_proc, int num_proc);

double single_spin_energy(int index, const upcxx::global_ptr<int>& grid, int row_size, double J, double B);

double energy(const upcxx::global_ptr<int>& grid, double J, double B, int row_size);

int* generateSpins(int rows_per_proc, int row_size, int rank);

void flipSpin(upcxx::global_ptr<int>& grid, int idx);

double avgMagnetism(const upcxx::global_ptr<int>& spinArray, int spinArraySize);

void saveGrid(const upcxx::global_ptr<int>& grid, int row_size, std::string folderName);

void saveMag(double mg, std::string folderName);

void saveEnergy(double energy, std::string folderName);

std::string createFolderWithTimestampName(int rep);

void saveParametersToFile(int netSize, double J, double B, long long iters, long long repeat);

void readParametersFromFile(int& netSize, double& J, double& B, long long& iters, long long& repeat);