#include "nvector.hpp"
#include "nvector_imp.hpp"
#include "matrix.hpp"
#include "matrix_imp.hpp"
#include "qr_decomposition.hpp"
#include <iostream>

int main(){
    // below is for test purpose

    // check inv, adjoint
    matrix<double> a(4, 4, { { 5, -2, 2, 7 }, { 1, 0, 0, 3 }, { -3, 1, 5, 0 },{ 3, -1, -9, 4 }});
//    a.display();
//    std::cout << std::endl;
//    a.adjoint().display();
//    std::cout << std::endl;
//    a.inv().display();
//    std::cout << std::endl;

    // check qr decomposition

    matrix<double> basis(3, 3, {{2.92, 0.86, -1.15}, {0.86, 6.51, 3.32}, {-1.15, 3.32, 4.57}});
//    std::pair<matrix<double>, matrix<double>> qr = qr_decomposition(basis, std::multiplies<double>());
//    std::cout << "Q: \n";
//    qr.first.display();
//    std::cout << "R: \n";
//    qr.second.display();

    // check eigen
    auto eig = eigen(basis, std::multiplies<double>());
    eig.first = eig.first.T();
    std::cout << "Eigenvectors: \n";
    eig.first.display();
    std::cout << "Eigenvalues: ";
    eig.second.display();
    std::cout << "Testing: \n";
    (basis * eig.first.mat[2]).display();
    (eig.first.mat[2] * eig.second[2]).display();
}
