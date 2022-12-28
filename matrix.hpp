#ifndef __matrix_hpp
#define __matrix_hpp

#include "nvector.hpp"
#include <vector>

#define MAXN 100

template<class type>
    class matrix{
        public:
            int N, M;
            std::vector<nvector<type>> mat;
            matrix(int Ni = 0, int Mi = 0, std::vector<nvector<type>> mati = {}){
                N = Ni;
                M = Mi;
                if (mati.size())
                    mat = mati;
                else mat.resize(N, nvector<type>(M, 0));
            }
            void display();
            matrix operator+(matrix b);
            matrix operator-(matrix b);
            matrix operator*(matrix b);
            nvector<type> operator*(nvector<type> b);
            matrix cofactor(int p, int q);
            type det();
            type trace();
            matrix adjoint();
            matrix inv();
            matrix pinv();
            matrix T();
            nvector<type> diag();
            bool is_upper_triag(double tol = 1e-9);
    };

#endif // __matrix_hpp
