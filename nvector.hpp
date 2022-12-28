#ifndef __nvector_hpp
#define __nvector_hpp

#include <vector>

template<class type>
    class nvector : private std::vector<type>{
        public:
            using std::vector<type>::vector;
            using std::vector<type>::push_back;
            using std::vector<type>::pop_back;
            using std::vector<type>::size;
            using std::vector<type>::resize;
            using std::vector<type>::begin;
            using std::vector<type>::end;
            using std::vector<type>::operator[];
            void display();
            double norm();
            double mean();
            nvector<type> operator+(nvector<type> b);
            nvector<type> operator-(nvector<type> b);
            nvector<double> operator*(double b);
            nvector<double> operator/(double b);
    };

#endif // __nvector_hpp
