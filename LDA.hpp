#include "matrix_imp.hpp"
#include "nvector_imp.hpp"
#include "qr_decomposition.hpp"

template <class type>
    matrix<type> vecxvec(nvector<type> a, nvector<type> b){
        matrix<type> res(a.size(), b.size());
        for (int i = 0; i < a.size(); i++)
            for (int j = 0; j < a.size(); j++)
                res.mat[i][j] = a[i] * b[j];
        return res;
    }

template <class type>
    matrix<double> Sw(matrix<type> data, nvector<int> class_assign){
        assert(data.N == class_assign.size());
        int num_class = *max_element(class_assign.begin(), class_assign.end());
        nvector<nvector<double>> mean(num_class + 1, nvector<double>(data.M));
        nvector<int> cnt(num_class + 1);
        for (int i = 0; i < class_assign.size(); i++){
            cnt[class_assign[i]]++;
            mean[class_assign[i]] = mean[class_assign[i]] + data.mat[i];
        }
        for (int i = 1; i <= num_class; i++)
            mean[i] = mean[i] / cnt[i];
        matrix<double> res(data.M, data.M);
        for (int i = 0; i < class_assign.size(); i++)
            res = res + vecxvec(data.mat[i] - mean[class_assign[i]], data.mat[i] - mean[class_assign[i]]);
        return res;
    }
template <class type>
    matrix<double> Sb(matrix<type> data, nvector<int> class_assign){
        assert(data.N == class_assign.size());
        int num_class = *max_element(class_assign.begin(), class_assign.end());
        nvector<nvector<double>> mean(num_class + 1, nvector<double>(data.M));
        nvector<int> cnt(num_class + 1);
        nvector<double> all(data.M);
        for (int i = 0; i < class_assign.size(); i++){
            cnt[class_assign[i]]++;
            mean[class_assign[i]] = mean[class_assign[i]] + data.mat[i];
            all = all + data.mat[i];
        }
        for (int i = 1; i <= num_class; i++)
            mean[i] = mean[i] / cnt[i];
        all = all / data.N;
        matrix<double> res(data.M, data.M);
        for (int i = 1; i <= num_class; i++)
            res = res + vecxvec(mean[i] - all, mean[i] - all) * data.N;
        return res;

    }
template <class type>
    matrix<type> LDA_decomposition(matrix<type> data, nvector<int> class_assign, int num_class){
        assert(data.N == class_assign.size());
        assert(num_class > 0);
        assert(num_class < data.M);
        assert(num_class < *max_element(class_assign.begin(), class_assign.end()));
        auto U = Sw(data, class_assign).inv() * Sb(data, class_assign);
        auto eigvec = eigen(U, std::multiplies<double>()).first;
        auto eigval = eigen(U, std::multiplies<double>()).second;
        eigvec = eigvec.T();
        eigvec.N = num_class;
        while (eigvec.mat.size() > eigvec.N)
            eigvec.mat.pop_back();
        eigvec = eigvec.T();
        return data * eigvec;
    }
