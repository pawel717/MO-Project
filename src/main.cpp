#include <iostream>
#include <fstream>
#include <cmath>
#include "calerf.h"
#include "Vector.h"
#include "Matrix.h"
#include "TridiagonalMatrix.h""
using namespace std;

void errorsRespectToH(Matrix& U_thomas, Matrix& U_gs, Matrix& U_analytic);
void errorsRespectToT(Matrix& U_thomas, Matrix& U_gs, Matrix& U_analytic);
void compareThomasGS(Matrix& U_thomas, Matrix& U_gs, Matrix& U_analytic);

/* Parametry wynikajace z postawionego problemu */
double D=1;
double t0=0;                // Poczatek przedzialu czasowego t
double t_max=2;             // Koniec przedzialu czasowego t
double x0=0;                // Poczatek przedzialu przestrzennego x
double x_max=6*sqrt(2);     // Koniec przedzialu przestrzennego x
double alfa=1;              // **** //
double beta=0;
double gamma=-1;            // Z warunków brzegowych
double fi=1;
double psi=0;
double eta=0;               //
double h;                   // Krok przestrzenny
double delta_t=h*h;

int main()
{
    Matrix U_thomas;        // Macierz z wynikami algorytmem Thomasa
    Matrix U_gs;            // Macierz z wynikami algorytmem Gaussa Seidela
    Matrix U_analytic;      // Macierz z wynikami analityczne

    errorsRespectToH(U_thomas, U_gs, U_analytic);
    errorsRespectToT(U_thomas, U_gs, U_analytic);
    compareThomasGS(U_thomas, U_gs, U_analytic);

    return 0;
}

/* Zwraca rozwiazanie analityczne */
double analytic(double x, double t)
{
    return ERFCL(x/(2.0*sqrt(D*t)));
}

/* Funkcja zapisujaca linie z danymi do pliku */
void save_line(fstream &file, double *data, int n)
{
    for(int i=0;i<n;i++)
    {
        file.width(15);
        file<<*(data+i);
    }
    file<<endl;
}

/* Funkcja zapisujaca linie z tytulem do pliku */
void save_title(fstream &file, string *data, int n)
{
    for(int i=0;i<n;i++)
    {
        file.width(15);
        file<<*(data+i);
    }
    file<<endl;
}

/* Zwraca wspolczynnik lambda z metody Laasonen */
double lambda(double delta_t, double h)
{
    return (D*delta_t)/(h*h);
}

/* Tworzenie macierzy trójdiagonalnej A do metody Laasonen */
void generateMatrixA(TridiagonalMatrix& A, double h, double delta_t, long long int n)
{
    for(int i=1;i<n-1;i++)
    {
        A(i,i-1)=lambda(delta_t, h);
        A(i,i)=-1*(1+2*lambda(delta_t,h));
        A(i,i+1)=lambda(delta_t,h);
    }
}

/* Ustawia elementy macierzy zgodnie z warunkami brzegowymi */
void boundaryCondition(TridiagonalMatrix& A, Matrix& U, Vector& W, double h, double delta_t)
{
    long long int m=(t_max-t0)/delta_t;     // Ilosc poziomów czasowych
    long long int n=(x_max-x0)/h;           // Ilosc punktow przestrznnych
    A(0,0)=alfa-beta/h;
    A(0,1)=beta/h;
    A(n-1,n-2)=-psi/h;
    A(n-1,n-1)=fi+psi/h;
    W(0)=-1*gamma;
    W(n-1)=-1*eta;
    for(int i=0;i<m;i++)
    {
        U(i,0)=1;
        U(i,n-1)=0;
    }
}

/* Ustawia elementy macierzy zgodnie z warunkami poczatkowymi */
void initialCondition(Matrix& U, long long int n)
{
    for(int i=0;i<n;i++) U(0,i)=0;
}

/* Tworzy równanie macierzowe zgodnie z algorytmem Laasonen AU=W */
void laasonen(TridiagonalMatrix& A, Matrix& U, Vector& W, double h, double delta_t)
{
    long long int m=(t_max-t0)/delta_t;     // Ilosc poziomów czasowych
    long long int n=(x_max-x0)/h;           // Ilosc punktow przestrznnych
    A=TridiagonalMatrix(n);
    U=Matrix(n,m);
    W=Vector(n);
    boundaryCondition(A,U,W,h,delta_t);  // warunki brzegowe
    initialCondition(U,n);                   // warunek poczatkowy
    generateMatrixA(A,h,delta_t,n);          //generacja macierzy A
}

/* Rozwiazuje rowanie rozniczkowe za pomoca rownan macierzowych
   z macierza trojdiagonalna algorytmem Thomasa AU_current=W
   i zapisuje wyniki w macierzy U */
void solveWithThomas(TridiagonalMatrix& A, Matrix& U, Vector& W)
{
    Vector U_current=Vector(U.get_size_n()); // Wektor na wyliczone rozwiazanie na poziomie tk+1
    A.thomas();                              // Dekompozycja macierzy A - wystarczy ją rozłożyć raz
    for(int k=0;k<U.get_size_m()-1;k++)
    {
        A.solve(W,U_current);               // Rozwiazanie rownania macierzowego AU=W
        U.setRow(U_current,k+1);            // Zapisanie wyliczonego wyniku na poziomie tk+1
        for(int i=1;i<U.get_size_n()-1;i++) W(i)=-U_current(i); // W=-Utk
    }
}

/* Rozwiazuje rowanie rozniczkowe za pomoca rownan macierzowych
   z macierza trojdiagonalna algorytmem Gaussa-Seidela Ax=b
   i zapisuje wyniki  w macierzy U */
void solveWithGaussSeidel(TridiagonalMatrix& A, Matrix& U, Vector& b) //(D+L)xn1=-Uxn +b
{
    double epsylon=1.0e-5;      // Dopuszczalny blad residuum
    double delta=1.0e-15;       // Dopuszczalny estymator bledu
    double estymator,residuum;
    Vector xn(U.get_size_n());  // Poprzednie przyblizenie
    Vector xn1(U.get_size_n()); // Aktualne przyblizenie
    int N=50;                   // Ilosc iteracji
    int n=A.get_size_n();

    for(int k=0;k<U.get_size_m()-1;k++)
    {
        for(int k=0;k<N;k++)   // N ilosc iteracji - pierwsze kryterium zatrzymania
        {
            for(int i=0;i<n;i++)
            {
                double sum1=0, sum2=0;
                for(int j=0;j<i;j++)
                {
                    sum1+=A(i,j)*xn1(j);
                }
                for(int j=i+1;j<n;j++)
                {
                    sum2+=A(i,j)*xn(j);
                }
                xn1(i)=(b(i)-sum1-sum2)/A(i,i);                // xn1 kolejne przyblizenie
            }
            residuum=(multiplyMatrixVector(A,xn1)-b).norm();   // Residuum - Ax-b=0
            estymator=(xn1-xn).norm();                         // Estymator - roznica miedzy przyblizeniami
            xn=xn1;
            if(estymator<delta || residuum<epsylon) break;     // Sprawdzenie kryteriow zatrzymania
        }
        for(int i=1;i<U.get_size_n()-1;i++) b(i)=-xn(i);       // W=-Utk
        U.setRow(xn,k+1);
    }
}

/* Rozwiazuje rowanie rozniczkowe za pomoca wzoru analitycznego
   i zapisuje wyniki  w macierzy U */
void solveWithAnalytic(Matrix& U, double delta_t, double h)
{
    long long int n=(x_max-x0)/h;          // Ilosc punktow przestrznnych
    long long int m=(t_max-t0)/delta_t;    // Ilosc poziomów czasowych
    U=Matrix(n,m);
    double t=t0;

    for(int j=0;j<m;j++)
    {
        double x=x0;
        for(int i=0;i<n;i++)
        {
            U(j,i)=analytic(x,t);
            x+=h;
        }
        t+=delta_t;
    }
}

/* Wyznacza i zapisuje do pliku bledy rozwiazan metodami
   Thomasa i Gaussa-Seidela w zaleznosci od kroku przestrzennego h */
void errorsRespectToH(Matrix &U_thomas, Matrix &U_gs, Matrix& U_analytic)
{
    fstream wyniki;
    wyniki.open("wyniki_bledy_od_h.txt", ios::out|ios::trunc);
    if(wyniki.good()==false) perror("Nie udalo sie otworzyc pliku");
    TridiagonalMatrix A;
    Vector W_thomas;
    Vector W_gs;
    h=1.0;

    for(int i=0;i<7;i++,h/=2)
    {
        delta_t=h*h;
        long long int n=(x_max-x0)/h;               // Ilosc punktow przestrznnych
        long long int m=(t_max-t0)/delta_t;         // Ilosc poziomów czasowych
        laasonen(A,U_thomas,W_thomas,h,delta_t);    // Stworzenie rowania macierzowego Laasonen AU=W
        W_gs=W_thomas;
        U_gs=U_thomas;
        solveWithGaussSeidel(A,U_gs,W_gs);          // Wyznaczenie rozwiazania metoda Gaussa-Seidela
        solveWithThomas(A,U_thomas,W_thomas);       // Wyznaczenie rozwiazania metoda Thomasa
        solveWithAnalytic(U_analytic,delta_t,h);    // Wyznaczenie rozwiazania analityczne
        double error_thomas=(U_analytic.getRow(m-1)-U_thomas.getRow(m-1)).norm(); // Blad- norma maximum
        double error_gs=(U_analytic.getRow(m-1)-U_gs.getRow(m-1)).norm();         // Blad- norma maximum
        double data[3]={log10(h),log10(error_thomas),log10(error_gs)};
        save_line(wyniki, data, 3);
        cout<<h<<" " <<delta_t<<endl;
    }

    h*=2;
    wyniki.close();
}

/* Wyznacza i zapisuje do pliku bledy rozwiazan metodami
   Thomasa i Gaussa-Seidela w zaleznosci od poziomu czasowego t */
void errorsRespectToT(Matrix& U_thomas, Matrix& U_gs, Matrix& U_analytic)
{
    fstream wyniki;
    wyniki.open("wyniki_bledy_od_t.txt", ios::out|ios::trunc);
    if(wyniki.good()==false) perror("Nie udalo sie otworzyc pliku");
    double t=t0;
    long long int n=(x_max-x0)/h;                                             // Ilosc poziomow czasowych

    for(int i=0;i<U_thomas.get_size_m();i++,t+=delta_t)
    {
        double error_thomas=(U_analytic.getRow(i)-U_thomas.getRow(i)).norm(); // Blad- norma maximum
        double error_gs=(U_analytic.getRow(i)-U_gs.getRow(i)).norm();         // Blad- norma maximum
        double data[3]={t,error_thomas,error_gs};
        save_line(wyniki, data, 3);
    }
    wyniki.close();
}

/* Zapisuje do pliku porownanie wynikow otrzymanych metodami Thomasa,
   Gaussa-Seidela i analityczne na poziomach czasowych 0.2, 0.9, 1.6 */
void compareThomasGS(Matrix& U_thomas, Matrix& U_gs, Matrix& U_analytic)
{
    fstream wyniki;
    wyniki.open("wyniki_analytic.txt", ios::out|ios::trunc);
    if(wyniki.good()==false) perror("Nie udalo sie otworzyc pliku");
    fstream wyniki2;
    wyniki2.open("wyniki_thomas.txt", ios::out|ios::trunc);
    if(wyniki2.good()==false) perror("Nie udalo sie otworzyc pliku");
    fstream wyniki3;
    wyniki3.open("wyniki_GS.txt", ios::out|ios::trunc);
    if(wyniki3.good()==false) perror("Nie udalo sie otworzyc pliku");
    long long int n=(x_max-x0)/h;   // Ilosc punktow przestrzennych
    double t[3]={0.2, 0.9, 1.6};    // Poziomy czasowe na ktorych porownujemy wyniki
    cout<<h<<endl;
    for(int i=0;i<3;i++)
    {
        wyniki<<t[i]<<endl;
        wyniki2<<t[i]<<endl;
        wyniki3<<t[i]<<endl;
        double x=x0;
        for(int j=0;j<n;j++)
        {
            double data[2]={x, U_analytic(t[i]/delta_t,j)};
            save_line(wyniki,data,2);
            double data2[2]={x, U_thomas(t[i]/delta_t,j)};
            save_line(wyniki2,data2,2);
            double data3[2]={x, U_gs(t[i]/delta_t,j)};
            save_line(wyniki3,data3,2);
            x+=h;
        }
    }
    wyniki.close();
    wyniki2.close();
    wyniki3.close();
}
