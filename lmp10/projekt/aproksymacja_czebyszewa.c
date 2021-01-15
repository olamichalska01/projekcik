#include <stdio.h>

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
		return (2*cze_n(st-1,x) + 2*x*d_czen_jeden(st-1,x)-d_czen_jeden(st-2,x);
	}	
}

double d_czen_dwa(int st, double x)
{
	return 0;
}

double d_czen_trzy(int st, double x)
{
	return 0;
}
