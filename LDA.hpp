#include "matrix_imp.hpp"
#include "nvector_imp.hpp"
#include "qr_decomposition.hpp"


template <class type>
    matrix<type> LDA_decomposition(matrix<type> data, nvector<int> class_assign, int num_class){
        assert(num_class < data.M);
        assert(num_class < *max_element(class_assign.begin(), class_assign.end()));
        return data;
    }
