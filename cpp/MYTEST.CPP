#include <iostream.h>
 #include <stdio.h>

#include <math.h>
void main()
{
	char *s;
	long double x;

	x = acos(30.0);

	printf("%3d\n", x );

}
 /*
long double CSC(long double XVal)
{
  static long double Infinity = 999999.9, Small = 0.000001;
  long double Temp;

  Temp = sin(XVal);

  if (fabs(Temp) < Small)
	 return Infinity;
  else
	 return (1.0 / Temp);
}  /* CSC */

