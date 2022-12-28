#ifndef __nvector_imp_hpp
#define __nvector_imp_hpp

#include "math.h"
#include "nvector.hpp"
#include <cassert>
#include <iostream>

template <class type>
    void nvector<type>::display(){
        for (auto i : *this)
            std::cout << i << " ";
        std::cout << std::endl;
    }

template <class type>
    double nvector<type>::norm(){
        double res = 0;
        for (auto i : *this)
            res += i * i;
        return std::sqrt(res);
    }
template <class type>
    double nvector<type>::mean(){
        double res = 0;
        for (auto i : *this)
            res += i;
        return res / this->size();
    }

template <class type>
    nvector<type> nvector<type>::operator+(nvector<type> b){
        assert(this->size() == b.size());
        nvector<type> res(b.size());
        for (int i = 0; i < b.size(); i++)
            res[i] = b[i] + (*this)[i];
        return res;
    }
template <class type>
    nvector<type> nvector<type>::operator-(nvector<type> b){
        assert(this->size() == b.size());
        nvector<type> res(b.size());
        for (int i = 0; i < b.size(); i++)
            res[i] = -b[i] + (*this)[i];
        return res;
    }
template <class type>
    nvector<double> nvector<type>::operator*(double b){
        nvector<double> res(this->size());
        for (int i = 0; i < this->size(); i++)
            res[i] = (*this)[i] * b;
        return res;
    }
template <class type>
    nvector<double> nvector<type>::operator/(double b){
        nvector<double> res(this->size());
        for (int i = 0; i < this->size(); i++)
            res[i] = (*this)[i] / b;
        return res;
    }
#endif // __nvector_imp_hpp
