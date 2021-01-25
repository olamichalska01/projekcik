#include <stdio.h>
#include <math.h>
#include "makespl.h"
#include "piv_ge_solver.h"
#include <stdlib.h>
#include <float.h>

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

        if(st == 0 || st == 1) return 0;
        else return -st * licznik / mianownik;
}

double d_czen_cztery(int k, double x)
{
        const double a = sqrt(x*x - 1);
        const double b = 6 * pow( k, 2 );
        const double c = x * x - 1;

        if(k == 0) return 0;

        return -(k * (pow((a+x), k) * (3*k*pow(c, 3/2)+(b+6)*x*x*x+a*((-k*k*k-14*k)*x*x+k*k*k-k)+(9-b)*x)+pow((x-a),k)*(3*k*pow(c, 3/2)+(-b-6)*x*x*x+a*((-k*k*k-14*k)*x*x+k*k*k-k)+(b-9)*x))) / (2*pow(c, 7/2));

}

void make_spl(points_t * pts, spline_t * spl)
{

        matrix_t       *eqs= NULL;
        double         *x = pts->x;
        double         *y = pts->y;
        double          a = x[0];
        double          b = x[pts->n - 1];
        int             i, j, k;
        int             nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  
	char *nbEnv= getenv( "APPROX_BASE_SIZE" );

        if( nbEnv != NULL && atoi( nbEnv ) > 0 )
                nb = atoi( nbEnv );

        eqs = make_matrix(nb, nb + 1);


  
	for (j = 0; j < nb; j++) 
	{
                for (i = 0; i < nb; i++)
                        for (k = 0; k < pts->n; k++)
                                add_to_entry_matrix(eqs, j, i, cze_n(i, x[k]) * cze_n(j, x[k]));

                for (k = 0; k < pts->n; k++)
                        add_to_entry_matrix(eqs, j, nb, y[k] * cze_n(j, x[k]));
        }
 
	if (piv_ge_solver(eqs)) 
	{
                spl->n = 0;
                return;
        }
 
	if (alloc_spl(spl, nb) == 0) 
	{
                for (i = 0; i < spl->n; i++) 
		{
                        double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
                        xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
                        spl->f[i] = 0;
                        spl->f1[i] = 0;
                        spl->f2[i] = 0;
                        spl->f3[i] = 0;
                        for (k = 0; k < nb; k++) 
			{
                                double          ck = get_entry_matrix(eqs, k, nb);
                                spl->f[i]  += ck * cze_n(k, xx);
                                spl->f1[i] += ck * d_czen_jeden(k, xx);
                                spl->f2[i] += ck * d_czen_dwa(k, xx);
                                spl->f3[i] += ck * d_czen_trzy(k, xx);
                        }
                }
        }
 }
