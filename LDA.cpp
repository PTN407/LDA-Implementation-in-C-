#include "LDA.hpp"
#include <iostream>

int main(){
    matrix<double> data(17, 2, {
        {2, 16},
        {4.46, 15.37},
        {6.48, 15.97},
        {4.3, 13.67},
        {2.44, 13.59},
        {1.96, 6.48},
        {3.17, 7.89},
        {4.87, 6.6},
        {3.01, 4.74},
        {3.17, 6.6},
        {0.95, 5.02},
        {-6.57, 11.49},
        {-7.22, 9.14},
        {-5.92, 7.61},
        {-5.76, 9.47},
        {-6.53, 10.44},
        {-7.38, 7.81}
        });
    nvector<int> classes = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3};
    auto new_data = LDA_decomposition(data, classes, 1);
    new_data.display();
//    proj(nvector<double>({1, -1}), nvector<double>({2, 1}), std::multiplies<double>()).display();

}
