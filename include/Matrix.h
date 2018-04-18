#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <fstream>
#include "Vector.h"
class Matrix
{
    private:
        int size_n=0;
        int size_m=0;
        std::vector<std::vector<double> > matrix;
    public:
        Matrix();
        Matrix(int, int);
        virtual ~Matrix();
        int get_size_n();
        int get_size_m();
        int maxInColumn(int n, int m);
        void switchRows(int k, int l);
        void init(double*);
        void print();
        Vector getColumn(int m);
        Vector getRow(int n);
        void setColumn(Vector column, int m);
        void setRow(Vector row, int n);
        //operacje na macierzach
        double norm();
        double norm2();
        double determinant();
        double minor(int i, int j);
        Matrix transpose();
        Matrix complement();
        Matrix removeColumnRow(int k, int l);
        Matrix scalarProduct(double x);
        Matrix absolute();
        Matrix invert();
        void saveToFile(std::fstream &file);
        //funkcje zaprzyjaŸnione
        friend Matrix multiplyMatrix(Matrix &matrix1, Matrix &matrix2);
        friend Vector multiplyMatrixVector(Matrix &matrix, Vector &vect);
        //prze³adowania
        Matrix operator+(const Matrix &m);
        Matrix operator-(const Matrix &m);
        bool operator<(const double &x);
        double& operator()(const int &n, const int &m);
};

#endif // MATRIX_H
