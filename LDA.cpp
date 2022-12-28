#include "LDA.hpp"
#include <iostream>

int main(){
    matrix<double> data(4, 2, {{1, 2}, {1, 3}, {2, 4}, {2, 5}});
    nvector<int> classes = {1, 1, 2, 2};
    auto new_data = LDA_decomposition(data, classes, 1);
    new_data.display();
}
