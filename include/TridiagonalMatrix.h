#ifndef TRIDIAGONALMATRIX_H
#define TRIDIAGONALMATRIX_H
#include "Vector.h"
#include <vector>

class TridiagonalMatrix
{
    private:
        static double zero;
        int size_n=0;
        int size_m=3;
        std::vector<std::vector<double> > matrix;
    public:
        TridiagonalMatrix();
        TridiagonalMatrix(int);
        virtual ~TridiagonalMatrix();
        int get_size_n();
        int get_size_m();
        void init(double*);
        void init(double*, double*, double*);
        void print();
        //operacje na macierzach
        void thomas();
        void solve(Vector &b, Vector &X);
        //przeladowania
        double& operator()(const int &n, const int &m);
        //funkcje zaprzyjaznione
        friend Vector multiplyMatrixVector(TridiagonalMatrix &matrix, Vector &vect);
};

#endif // TRIDIAGONALMATRIX_H
