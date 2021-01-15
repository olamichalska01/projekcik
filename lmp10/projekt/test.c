#include <stdio.h>
#include <math.h>

double cze_n(int st, double x)
{
        if(st == 0) return 1;
        else if( st == 1) return x;
        else  return 2*x*cze_n(st-1,x)-cze_n(st-2,x);
}

double d_czen_jeden(int st, double x)
{
        if(st == 0) return 0;
        else if( st == 1) return 1;
        else
        {
                return (2*cze_n(st-1,x) + 2*x*d_czen_jeden(st-1,x)-d_czen_jeden(st-2,x));
        }
}


double d2(int st, double x)
{
	return st/2*( pow( x + pow( x*x -1, 1/2) , st - 1 ) * ( 1 + x * pow(x*x-1, -1/2) ) + pow( x - pow( x*x -1, 1/2) , st - 1 ) * ( 1 - x * pow(x*x-1, -1/2) ) );
}

int main()
{

	int s = 3;
	double x = 5.4;

	printf("\t %lf \n\t %lf\n",d_czen_jeden(s,x),d2(s,x)); 

	return 0;
}
