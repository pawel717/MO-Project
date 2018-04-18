#include "C:/Users/Pawel1/Desktop/METODY OBLICZENIOWE/LABORATORIA/projekt/include/TridiagonalMatrix.h"
#include <cmath>
#include <iostream>
using namespace std;
double TridiagonalMatrix::zero=0;

TridiagonalMatrix::TridiagonalMatrix() {}

TridiagonalMatrix::TridiagonalMatrix(int size_n)
{
    this->size_n=size_n;
    this->matrix=vector<vector<double> >(this->size_m,vector<double>());
    this->matrix[0]=vector<double>(size_n-1,0);
    this->matrix[1]=vector<double>(size_n,0);
    this->matrix[2]=vector<double>(size_n-1,0);
    //ctor
}

TridiagonalMatrix::~TridiagonalMatrix()
{
    matrix[0].capacity();
    matrix[1].capacity();
    matrix[2].capacity();
    matrix.capacity();
}

void TridiagonalMatrix::init(double* elements)
{
    for(int i=0;i<this->size_m;i++)
    {
        for(int j=0;j<this->matrix[i].size();j++)
        {
            this->matrix[i][j]=*(elements++);
        }
    }
}

void TridiagonalMatrix::init(double* upper, double* diagonal, double* lower)
{
    int i=0;
    for(i=0;upper[i]!=NULL;i++)
    {
        this->matrix[0][i]=upper[i];
        this->matrix[1][i]=diagonal[i];
        this->matrix[2][i]=lower[i];
    }
    this->matrix[1][i]=diagonal[i];
}

// Gettery
int TridiagonalMatrix::get_size_n()
{
    return this->size_n;
}

int TridiagonalMatrix::get_size_m()
{
    return this->size_n;
}

//przeladowania operatorow
double& TridiagonalMatrix::operator()(const int &n, const int &m)
{
    if(n>this->size_n || n<0 || m>this->size_n || m<0) return zero;
    if(m==n+1) return this->matrix[0][m-1];
    if(m==n) return this->matrix[1][m];
    if(m+1==n) return this->matrix[2][m];
    return zero;
}
void TridiagonalMatrix::print()
{
    for(int i=0;i<this->size_n;i++)
    {
        for(int j=0;j<this->size_n;j++)
        {
            cout.width(9);
            cout<<(*this)(i,j)<<"     ";
        }
        cout<<endl;
    }
}

/* Algorytm thomasa - dekompozycja macierzy A */
void TridiagonalMatrix::thomas()
{
    for(int i=1;i<this->size_n;i++)
    {
        this->matrix[1][i]=this->matrix[1][i]-(this->matrix[2][i-1]/this->matrix[1][i-1]*this->matrix[0][i-1]);
    }
}

/* Algorytm Thomasa - rozwiazanie ukladu poprzez podstawienie wsteczne */
void TridiagonalMatrix::solve(Vector &b, Vector &X)
{
    for(int i=1;i<this->size_n;i++)
    {
        b(i)=b(i)-(this->matrix[2][i-1]/this->matrix[1][i-1]*b(i-1)); // Podstawienie w przód
    }
    X(size_n-1)=b(size_n-1)/this->matrix[1][size_n-1];
    for(int i=size_n-2;i>=0;i--)
    {
        X(i)=(b(i)-(this->matrix[0][i]*X(i+1)))/this->matrix[1][i];   // Podstawienie wstecz
    }
}

//funkcje zaprzyjaznione
Vector multiplyMatrixVector(TridiagonalMatrix &matrix, Vector &vect)
{
    int size_k=vect.get_size();
    int size_m=matrix.size_m;
    int size_n=matrix.size_n;
    Vector newVect(size_m);
        for(int i=0;i<size_m;i++)
        {
            double sum=0;
            for(int j=0;j<size_n;j++)
            {
                sum+=matrix(i,j)*vect(j);
            }
            newVect(i)=sum;
        }
    return newVect;
}
