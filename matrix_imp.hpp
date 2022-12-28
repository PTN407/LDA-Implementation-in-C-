#include "matrix.hpp"
#include "nvector.hpp"
#include "nvector_imp.hpp"
#include <cassert>
#include <functional>
#include <iostream>

template <class type>
    void matrix<type>::display(){
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->M; j++)
                std::cout << this->mat[i][j] << "\t\n"[j == this->M - 1];
    }
template <class type>
    matrix<type> matrix<type>::operator+(matrix<type> b){
        assert((this->N == b.N) && (this->M == b.M));
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->M; j++)
                b.mat[i][j] += this->mat[i][j];
        return b;
    };
template <class type>
    matrix<type> matrix<type>::operator-(matrix<type> b){
        assert((this->N == b.N) && (this->M == b.M));
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->M; j++)
                b.mat[i][j] = this->mat[i][j] -  b.mat[i][j];
        return b;
    };
template <class type>
    matrix<type> matrix<type>::operator*(matrix<type> b){
        assert(this->M == b.N);
        matrix<type> res(this->N, b.M);
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->M; j++)
                for (int k = 0; k < b.M; k++)
                    res.mat[i][k] += this->mat[i][j] *  b.mat[j][k];
        return res;
    };
template <class type>
    nvector<type> matrix<type>::operator*(nvector<type> b){
        assert(this->M == b.size());
        nvector<type> res(this->N);
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->M; j++)
                    res[i] += this->mat[i][j] * b[j];
        return res;
    };
template <class type>
    matrix<type> matrix<type>::T(){
        matrix<type> res(this->M, this->N);
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->M; j++)
                res.mat[j][i] = this->mat[i][j];
        return res;
    };
template <class type>
    matrix<type> matrix<type>::cofactor(int p, int q){
        matrix<type> res(this->N - 1, this->M - 1);
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->M; j++)
                if ((i != p && j != q))
                    res.mat[i - (i >= p)][j - (j >= q)] = this->mat[i][j];
        return res;
    }
template <class type>
    type matrix<type>::det(){
        assert(this->N == this->M);
        if (this->N == 1)
            return this->mat[0][0];
        matrix<type> temp;
        int sign = 1;
        type res = 0;
        for (int f = 0; f < this->N; f++){
            temp = this->cofactor(0, f);
            res += sign * this->mat[0][f] * temp.det();
            sign = -sign;
        }
        return res;
    }
template <class type>
    type matrix<type>::trace(){
        assert(this->N == this->M);
        double res = 0;
        for (int f = 0; f < this->N; f++)
            res += *this[f][f];
        return res;
    }

template <class type>
    nvector<type> matrix<type>::diag(){
        assert(this->N == this->M);
        nvector<type> res;
        for (int f = 0; f < this->N; f++)
            res.push_back(this->mat[f][f]);
        return res;
    }
template <class type>
    bool matrix<type>::is_upper_triag(double tol){
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < i; j++)
                if (abs(this->mat[i][j]) > tol) return 0;
        return 1;
    }
template <class type>
    matrix<type> matrix<type>::adjoint(){
        assert(this->N == this->M);
        matrix<type> res(this->N, this->N);
        int sign = 1;
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->N; j++){
                matrix<type> temp = this->cofactor(i, j);
                sign = ((i + j) % 2 == 0 ? 1 : -1);
                res.mat[j][i] = sign * temp.det();
            }
        return res;
    }
template <class type>
    matrix<type> matrix<type>::inv(){
        assert(this->N == this->M);
        type determinant = this->det();
        if (determinant == 0)
            return this->pinv();
        matrix<type> adj = this->adjoint();
        matrix<type> res(this->N, this->N);
        for (int i = 0; i < this->N; i++)
            for (int j = 0; j < this->N; j++)
                res.mat[i][j] = adj.mat[i][j] / determinant;
        return res;
    }
template <class type>
    matrix<type> matrix<type>::pinv(){
        matrix<type> a = (this->T() * (*this));
        matrix<type> b = this->T();
        return a.inv() * b;
    }

