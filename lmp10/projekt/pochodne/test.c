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
	if(st == 0) return 0;	
	else if(st == 1) return 1;
	else return 0.5 * st * ( pow( x + pow(x * x - 1, 0.5) , st - 1 ) * ( 1 + x * pow(x * x - 1, -0.5) ) + pow( x - pow(x * x - 1, 0.5) , st - 1 ) * ( 1 - x * pow(x * x - 1, -0.5) ) );
}

double d_czen_dwa(int st, double x) // jak ktos idzie do piekla to tam liczenie kolejnych pochodnych jest jedna z tortur
{

	if(st == 0 || st == 1) return 0;
       	

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


int main()
{

	int s;
	double x = 5.4;


	s = 0;
        printf("\ndla stopnia 0:\n\npierwsza:    %lf \ndruga:       %lf\ntrzecia:     %lf\nczwarta:     %lf\nczebyszew:   %lf\n", d_czen_jeden(s,x), d_czen_dwa(s,x), d_czen_trzy(s,x), d_czen_cztery(s,x), cze_n(s,x) );
	
	s = 1;
	printf("\ndla stopnia 1:\n\npierwsza:    %lf \ndruga:       %lf\ntrzecia:     %lf\nczwarta:       %lf\nczebyszew:   %lf\n", d_czen_jeden(s,x), d_czen_dwa(s,x), d_czen_trzy(s,x), d_czen_cztery(s,x), cze_n(s,x) ); 

	s = 3;
	printf("\ndla innego stopnia:\n\npierwsza:    %lf \ndruga:       %lf\ntrzecia:     %lf\nczwarta:      %lf\nczebyszew:   %lf\n\n\n", d_czen_jeden(s,x), d_czen_dwa(s,x), d_czen_trzy(s,x), d_czen_cztery(s,x), cze_n(s,x) );


	return 0;
}
