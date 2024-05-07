using System;
using static System.Math;

public class integration{
	public static (double,int) integrate
	(Func<double,double> f, double a, double b,
	double δ=0.00001, double ε=0.00001, double f2=Double.NaN, double f3=Double.NaN)
	{
		int evals = 0;

		double h=b-a;
		if(Double.IsNaN(f2)){ f2=f(a+2*h/6); f3=f(a+4*h/6); evals+=2;} // first call, no points to reuse
		double f1=f(a+h/6), f4=f(a+5*h/6); evals+=2;
		double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // higher order rule
		double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule
		double err = Math.Abs(Q-q);
		if (err <= δ+ε*Math.Abs(Q)) return (Q,evals);
		else {
			(double int1, int evals1) = integrate(f,a,(a+b)/2,δ/Math.Sqrt(2),ε,f1,f2);
			(double int2, int evals2) = integrate(f,(a+b)/2,b,δ/Math.Sqrt(2),ε,f3,f4);
			return (int1+int2,evals+evals1+evals2);
		}
	}
	public static (double,int) ccIntegrate
	(Func<double,double> f, double a, double b,
	double δ=0.001, double ε=0.001)
	{
		Func<double,double> newfunc = x => f( (a+b)/2+(b-a)/2*Math.Cos(x) )*Math.Sin(x)*(b-a)/2;
		return integrate(newfunc,0,Math.PI);
	}
}
