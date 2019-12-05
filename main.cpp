#include <iostream>
#include "../Lib/Math.h"
using namespace std;

inline double Diff_u1 (const Matrix<double>& X)
{
    assert (X.Get_Size() == 3);
    double x = *(X.Get_pointer()+2);
    return (-*X.Get_pointer()**(X.Get_pointer()+1) + x > 1e-10 ? sin(x)/x : 1);
}
inline double Diff_u2 (const Matrix<double>& X)
{
    assert (X.Get_Size() == 3);
    double x = *(X.Get_pointer()+2);
    return (-Pow(*(X.Get_pointer()+1), 2) + 3.125*x/(1+x*x));
}
int main()
{
    Array_of_Functions2 A (2);
    A[0] = Diff_u1;
    A[1] = Diff_u2;
    Matrix<double> u0 (1, 2);
    u0[0] = 0;
    u0[1] = -0.412;
    Matrix<double> Eps (1, 2);
    Eps[1] = Eps[0] = 1e-3;
    Solve_Differential_Equations::Explicit_Euler_method (A, 0, 1, u0, 10, Eps);
    return 0;
}
