#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

float f_float(float x)
{
    return sin(x) + cos(3*x);
}

float df_float(float x,float h)
{
    return (f_float(x+h)-f_float(x))/h;
}


double f_double(double x)
{
    return sin(x) + cos(3*x);
}

double df_double(double x,double h)
{
    return (f_double(x+h)-f_double(x))/h;
}

long double f_long_double(long double x)
{
    return sin(x) + cos(3*x);
}

long double df_long_double(long double x,long double h)
{
    return (f_long_double(x+h)-f_long_double(x))/h;
}

long double df_real_long_double(long double x)
{
    return cos(x) - 3*sin(3*x);
}

double df_real_double(double x)
{
    return cos(x) - 3*sin(3*x);
}

float df_real_float(float x)
{
    return cos(x) - 3*sin(3*x);
}




int main()
{
    cout << "====================" << endl;
    cout << "float" << endl;  
    cout << "====================" << endl;

    cout << setprecision(22) << fixed;

    float x = 1.0;
    for(int i=0;i<=40;i++)
    {
        float h = pow(2,-i);
        cout << "i = " << i  << " h = " << h << " 1 + h " <<  1+h << " df = " << df_float(x,h) << " real df = " << df_real_float(x) << endl;

    }
    cout << "====================" << endl;

    cout << "double" << endl;

    cout << "====================" << endl;

    double x_double = 1.0;
    for(int i=0;i<=40;i++)
    {
        double h = pow(2,-i);
        cout << "i = " << i  <<  " h = " << h << " 1 + h " <<  1+h << " df = " << df_double(x_double,h) << " real df = " << df_real_double(x_double) << endl;

    }

    cout << "====================" << endl;
    cout << "long double" << endl;
    cout << "====================" << endl;

    long double x_long_double = 1.0;
    for(int i=0;i<=40;i++)
    {
        long double h = pow(2,-i);
        cout << "i = " << i  << " h = " << h << " 1 + h " <<  1+h << " df = " << df_long_double(x_long_double,h) << " real df = " << df_real_long_double(x_long_double) << endl;

    }




    return 0;   
}