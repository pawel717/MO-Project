#include "Vector.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

//konstruktor
Vector::Vector(){}

Vector::Vector(int size_n)
{
    this->vect=vector<double>(size_n,0);
}

Vector::~Vector()
{
    vect.clear();
}

void Vector::print()
{
    for(int i=0;i<this->get_size();i++)
    {
        cout.width(15);
        cout<<vect[i]<<"     ";
    }
    cout<<endl;
}

void Vector::saveToFile(fstream &file)
{
    for(int i=0;i<this->get_size();i++)
    {
        file.width(15);
        file<<vect[i];
    }
    file<<endl;
}

void Vector::init(double* elements)
{
    for(int i=0;i<this->get_size();i++)
    {
        this->vect[i]=*(elements++);
    }
}

//settery/gettery
int Vector::get_size()
{
    return this->vect.size();
}

// operacje na macierzach

int Vector::maxInVector(int n)
{
    double maximum=abs(this->vect[n]);
    double maxi=n;
    for(int i=n+1;i<this->get_size();i++)
    {
        if(abs(this->vect[i]>maximum))
        {
            maximum=abs(this->vect[i]);
            maxi=i;
        }
    }
    return maxi;
}

double Vector::norm()
{
    double maximum=abs(this->vect[0]);
    for(int i=0;i<this->get_size();i++)
    {
        if(abs(this->vect[i])>maximum) maximum=abs(this->vect[i]);
    }
    return maximum;
}

Vector Vector::scalarProduct(double x)
{
    Vector newVector(this->get_size());
    for(int i=0;i<this->get_size();i++)
    {
        newVector.vect[i]=this->vect[i]*x;
    }
    return newVector;
}

Vector Vector::absolute()
{
    Vector newVector(this->get_size());
    for(int i=0;i<this->get_size();i++)
    {
        newVector.vect[i]=abs(this->vect[i]);
    }
    return newVector;
}

//przeladowania operatorow
double& Vector::operator()(const int &n)
{
    return this->vect[n];
}
Vector Vector::operator+(const Vector &v)
{
    Vector newVector(this->get_size());
    for(int i=0;i<this->get_size();i++)
    {
        newVector.vect[i]=this->vect[i]+v.vect[i];
    }
    return newVector;
}

Vector Vector::operator-(const Vector &v)
{
    Vector newVector(this->get_size());
    for(int i=0;i<this->get_size();i++)
    {
        newVector.vect[i]=this->vect[i]-v.vect[i];
    }
    return newVector;
}

bool Vector::operator<(const double &x)
{
    for(int i=0;i<this->get_size();i++)
    {
        if(this->vect[i]>=x) return false;
    }
    return true;
}
