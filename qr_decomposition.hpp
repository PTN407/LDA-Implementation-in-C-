#include "matrix.hpp"
#include "nvector.hpp"
#include "nvector_imp.hpp"
#include <cassert>
#include <functional>

template <class type, class prod>
    type inner_product(nvector<type> a, nvector<type> b, prod func){
        type res = 0;
        assert(a.size() == b.size());
        for (int i = 0; i < a.size(); i++)
            res += func(a[i], b[i]);
        return res;
    }
template <class type, class prod>
    nvector<double> proj(nvector<type> a, nvector<type> b, prod func){
        return a * inner_product(a, b, func) / inner_product(a, a, func);
    }
template <class type, class prod>
    matrix<type> orthonormalize(matrix<type> basis, prod inner_prod){
        matrix<type> res;
        res.mat.push_back(basis.mat.front());
        res.M = basis.mat.front().size();
        res.N = 1;
        for (int i = 1; i < basis.mat.size(); i++){
            nvector<type> projection(basis.M, 0);
            for (int j = 0; j < i; j++)
                projection = projection - proj(res.mat[j], basis.mat[i], inner_prod);
            res.mat.push_back(basis.mat[i] + projection);
            ++res.N;
        }
        for (int i = 0; i < basis.mat.size(); i++)
            res.mat[i] = res.mat[i] / res.mat[i].norm();
        return res;
    }

template <class type, class prod>
    std::pair<matrix<type>, matrix<type>> qr_decomposition(matrix<type> a, prod inner_prod){
        a = a.T();
        auto orthonormal_matrix = orthonormalize(a, inner_prod);
        auto q = orthonormal_matrix.T();
        matrix<type> r(q.M, q.N);
        for (int i = 0; i < r.N; i++)
            for (int j = 0; j <= i; j++)
                r.mat[i][j] = inner_product(orthonormal_matrix.mat[j], a.mat[i], inner_prod);
        return std::make_pair(q, r.T());
    }
template <class type, class prod>
    std::pair<matrix<double>, nvector<double>> eigen(matrix<type> a, prod inner_prod, int max_iter = 500, double tol = 1e-18){
        auto qr = qr_decomposition(a, inner_prod);
        matrix<double> previous(qr.first.N, qr.first.M);
        matrix<double> eigen_vector = qr.first;
        auto X = qr.second * qr.first;
        for (int i = 0; i < max_iter; i++){
            previous = qr.first;
            X = qr.second * qr.first;
            qr = qr_decomposition(X, inner_prod);
            eigen_vector = eigen_vector * qr.first;
            if (X.is_upper_triag(tol))
                break;
        }
        return std::make_pair(eigen_vector, X.diag());
    }
