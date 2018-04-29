#include <omp.h>
#include <chrono>
#include <iostream>
#include "eigen.hpp"
using namespace std;

int main() {
    DenseMatrix<double> a = {
        {1, 1, 1},
        {1, 2, 1},
        {1, 1, 2},
    };
    auto b = Eigenvalue<double>(a);
    DenseMatrix<double> v = b.getV();
    cout << v * Transpose(v);
    // DenseMatrix<double> a(DenseMatrix<double>::Identity(1000));
    // auto start = std::chrono::high_resolution_clock::now();
    // auto b = a + a;
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> diff = end - start;
    // std::cout << "Time to fill and iterate a vector of "
    //           << " ints : " << diff.count() << " s\n";
    return 0;
}