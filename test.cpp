#include <upcxx/upcxx.hpp>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]){
    upcxx::init();
    int rank = upcxx::rank_me();
    int num_procs = upcxx::rank_n();
    upcxx::global_ptr<int> myGlobalPtr = upcxx::new_array<int>(1);
    if(rank == 0){
        upcxx::rput(1, myGlobalPtr).wait();
    }
    upcxx::broadcast(&myGlobalPtr, 1, 0).wait();
    upcxx::global_ptr<int> sharedPtr = myGlobalPtr;  // Use the obtained value
    int value = upcxx::rget(sharedPtr).wait();                          // Read the shared data
    std::cout << value << std::endl;                                // Modify the shared data
    upcxx::finalize();
    return 0;
}