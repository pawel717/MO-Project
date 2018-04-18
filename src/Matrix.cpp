#include "Matrix.h"
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

//konstruktor
Matrix::Matrix(){}

Matrix::Matrix(int size_n, int size_m)
{
    this->size_n=size_n;
    this->size_m=size_m;
    this->matrix=vector<vector<double> >(this->size_n,vector<double>(this->size_m,0));
}

Matrix::~Matrix()
{
    //dtor
}

void Matrix::print()
{
    for(int i=0;i<this->size_m;i++)
    {
        for(int j=0;j<this->size_n;j++)
        {
            cout.width(9);
            cout<<this->matrix[j][i]<<"     ";
        }
        cout<<endl;
    }
}

void Matrix::init(double* elements)
{
    for(int i=0;i<this->size_m;i++)
    {
        for(int j=0;j<this->size_n;j++)
        {
            this->matrix[j][i]=*(elements++);
        }
    }
}

//settery/gettery
int Matrix::get_size_n()
{
    return this->size_n;
}

int Matrix::get_size_m()
{
    return this->size_m;
}

Vector Matrix::getRow(int n)
{
    Vector row(this->size_n);
    double tmp[this->size_n];
    for(int i=0;i<this->size_n;i++)
    {
        tmp[i]=this->matrix[i][n];
    }
    row.init(tmp);
    return row;
}

Vector Matrix::getColumn(int m)
{
    Vector column(this->size_m);
    double tmp[this->size_m];
    for(int i=0;i<this->size_m;i++)
    {
        tmp[i]=this->matrix[m][i];
    }
    column.init(tmp);
    return column;
}

void Matrix::setRow(Vector row, int n)
{
    for(int i=0;i<this->size_n;i++)
    {
        this->matrix[i][n]=row(i);
    }
}

void Matrix::setColumn(Vector column, int m)
{
    for(int i=1;i<this->size_m-1;i++)
    {
        this->matrix[m][i]=column(i);
    }
}

// operacje na macierzach

void Matrix::switchRows(int k, int l)
{
    for(int i=0;i<this->size_n;i++)
    {
        double temp=this->matrix[i][k];
        this->matrix[i][k]=this->matrix[i][l];
        this->matrix[i][l]=temp;
    }
}

int Matrix::maxInColumn(int n,int m)
{
    double maximum=abs(this->matrix[n][m]);
    double maxi=m;
    for(int i=m+1;i<this->size_m;i++)
    {
        if(abs(this->matrix[n][i])>maximum)
        {
            maximum=abs(this->matrix[n][i]);
            maxi=i;
        }
    }
    return maxi;
}

double Matrix::norm()
{
    double maximum=0;
    for(int i=0;i<this->size_n;i++)
    {
        double sum=0;
        for(int j=0;j<this->size_m;j++)
        {
            sum+=abs(this->matrix[i][j]);
        }
        if(sum>maximum) maximum=sum;
    }
    return maximum;
}
double Matrix::norm2()
{
    double maximum=0;
    for(int i=0;i<this->size_m;i++)
    {
        double sum=0;
        for(int j=0;j<this->size_n;j++)
        {
            sum+=abs(this->matrix[j][i]);
        }
        if(sum>maximum) maximum=sum;
    }
    return maximum;
}

double Matrix::determinant()
{
    if(this->size_n!=this->size_m)
    {
        cout<<"Nie mozna policzyc wyznacznika - macierz nie jest kwadratowa"<<endl;
        return NULL;
    }
    if(this->size_n<2)
    {
        if(this->size_n==2)
        {
            double a=this->matrix[0][0]*this->matrix[1][1];
            double b=this->matrix[0][1]*this->matrix[1][0];
            return a-b;
        }
        return this->matrix[0][0];
    }
    double w=0;
    for(int i=0;i<size_n;i++)
    {
        Matrix newMatrix=this->removeColumnRow(i,0);
        int sign;
        if((i+2)%2) sign=-1;
        else sign=1;
        double v=sign*this->matrix[i][0]*newMatrix.determinant();
        w+=v;
    }
    return w;
}

double Matrix::minor(int k, int l)
{
    Matrix minorMatrix(this->size_n-1,this->size_m-1);
    minorMatrix=this->removeColumnRow(k,l);
    return minorMatrix.determinant();
}

Matrix Matrix::transpose()
{
    Matrix newMatrix(this->size_m,this->size_n);
    for(int i=0;i<this->size_n;i++)
    {
        for(int j=0;j<this->size_m;j++)
        {
            newMatrix.matrix[j][i]=this->matrix[i][j];
        }
    }
    return newMatrix;
}

Matrix Matrix::removeColumnRow(int k, int l)
{
    Matrix newMatrix(this->size_n-1,this->size_m-1);
    for(int i=0, ii=0;i<this->size_n;i++)
    {
        if(i==k) continue;
        for(int j=0, jj=0;j<this->size_m;j++)
        {
            if(j==l) continue;
            newMatrix.matrix[ii][jj]=this->matrix[i][j];
            jj++;
        }
        ii++;
    }
    return newMatrix;
}

Matrix Matrix::complement()
{
    int w;
    Matrix newMatrix(this->size_n,this->size_m);
    for(int i=0;i<size_n;i++)
    {
        for(int j=0;j<size_m;j++)
        {
            if((i+j+2)%2) w=-1; // wspolczynnik
            else w=1;
            newMatrix.matrix[i][j]=w*minor(i,j);//*minor
        }
    }
    return newMatrix;
}

Matrix Matrix::scalarProduct(double x)
{
    Matrix newMatrix(this->size_n,this->size_m);
    for(int i=0;i<this->size_n;i++)
    {
        for(int j=0;j<this->size_m;j++)
        {
            newMatrix.matrix[i][j]=this->matrix[i][j]*x;
        }
    }
    return newMatrix;
}

Matrix Matrix::invert()
{
    double det=this->determinant();
    Matrix comp=this->complement();
    Matrix trans=comp.transpose();
    trans=trans.scalarProduct(1./det);
    return trans;
}

Matrix Matrix::absolute()
{
    Matrix newMatrix(this->size_n, this->size_m);
    for(int i=0;i<this->size_n;i++)
    {
        for(int j=0;j<this->size_m;j++)
        {
            newMatrix.matrix[i][j]=abs(this->matrix[i][j]);
        }
    }
    return newMatrix;
}

//przeladowania operatorow
double& Matrix::operator()(const int &n, const int &m)
{
    return this->matrix[m][n];
}

Matrix Matrix::operator+(const Matrix &m)
{
    Matrix newMatrix(this->size_n, this->size_m);
    for(int i=0;i<this->size_n;i++)
    {
        for(int j=0;j<this->size_m;j++)
        {
            newMatrix.matrix[i][j]=this->matrix[i][j]+m.matrix[i][j];
        }
    }
    return newMatrix;
}

Matrix Matrix::operator-(const Matrix &m)
{
    Matrix newMatrix(this->size_n, this->size_m);
    for(int i=0;i<this->size_n;i++)
    {
        for(int j=0;j<this->size_m;j++)
        {
            newMatrix.matrix[i][j]=this->matrix[i][j]-m.matrix[i][j];
        }
    }
    return newMatrix;
}

bool Matrix::operator<(const double &x)
{
    for(int i=0;i<this->size_n;i++)
    {
        for(int j=0;j<this->size_m;j++)
        {
            if(this->matrix[i][j]>=x) return false;
        }
    }
    return true;
}

void Matrix::saveToFile(fstream &file)
{
    for(int i=0;i<this->size_m;i+=1)
    {
        for(int j=0;j<this->size_n;j+=1)
        {
            file.width(15);
            file<<this->matrix[j][i];
        }
        file<<endl;
    }
}

//funkcje zaprzyjaznione
Matrix multiplyMatrix(Matrix &matrix1, Matrix &matrix2)
{
    int size_k=matrix2.size_n;
    int size_m=matrix1.size_m;
    int size_n=matrix1.size_n;
    Matrix newMatrix(size_k, size_m);
    for(int k=0;k<size_k;k++)
    {
        for(int i=0;i<size_m;i++)
        {
            double sum=0;
            for(int j=0;j<size_n;j++)
            {
                sum+=matrix1.matrix[j][i]*matrix2.matrix[k][j];
            }
            newMatrix.matrix[k][i]=sum;
        }
    }
    return newMatrix;
}
Vector multiplyMatrixVector(Matrix &matrix, Vector &vect)
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
                sum+=matrix.matrix[j][i]*vect(j);
            }
            newVect(i)=sum;
        }
    return newVect;
}
