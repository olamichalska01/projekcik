#include <stdio.h>
#include <math.h>



double cze_n(int st, double x)
{
        if(st == 0) return 1;
        else if( st == 1) return x;
        else  return 2 * x * cze_n(st-1, x) - cze_n(st-2, x);
}

double d_czen_jeden(int st, double x)
{
        return 0.5 * st * ( pow( x + pow(x * x - 1, 0.5) , st - 1 ) * ( 1 + x * pow(x * x - 1, -0.5) ) + pow( x - pow(x * x - 1, 0.5) , st - 1 ) * ( 1 - x * pow(x * x - 1, -0.5) ) );

}

double d_czen_dwa(int st, double x) // jak ktos idzie do piekla to tam liczenie kolejnych pochodnych jest jedna z tortur
{
        const double t = pow( x * x - 1, -0.5) * x;
        const double b = sqrt(x * x - 1);

        // podzielone na dwa wyrazenia, bo roznia sie tylko kilkoma znakami
        const double w1 = (st - 1) * pow( x + b, st - 2 ) * (1 + t) * (1 + t) + pow( x + b, st - 1 ) * ( -0.5 * pow( x * x - 1, -1.5 ) * 2 * x * x + pow( x * x - 1, -0.5 ) );
        const double w2 = (st - 1) * pow( x - b, st - 2 ) * (1 - t) * (1 - t) + pow( x - b, st - 1 ) * (  0.5 * pow( x * x - 1, -1.5 ) * 2 * x * x - pow( x * x - 1, -0.5 ) );

        return 0.5 * st * (w1 + w2);
}

double d_czen_trzy(int st, double x)
{
        const double a = sqrt(x*x - 1);

        const double licznik = pow( x - a, st) * ( 3 * st * x * a + (st*st + 2) * x * x - st * st + 1 ) + pow( x + a, st) * ( 3 * st * x * a + (-1 * st * st - 2)* x * x + st*st - 1 );

        const double mianownik = 2 * pow( x*x - 1, 2.5);


        return -st * licznik / mianownik;
}

