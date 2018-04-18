#ifndef VECTOR_H
#define VECTOR_H
#include <fstream>
#include <vector>

class Vector
{
private:
    std::vector<double> vect;
public:
    Vector();
    Vector(int);
    virtual ~Vector();
    int get_size();
    int maxInVector(int n);
    void init(double*);
    void print();
    void saveToFile(std::fstream &file);
    //operacje na macierzach
    double norm();
    Vector scalarProduct(double x);
    Vector absolute();
    //prze³adowania
    Vector operator+(const Vector &m);
    Vector operator-(const Vector &m);
    bool operator<(const double &x);
    double& operator()(const int &n);
};

#endif // VECTOR_H
