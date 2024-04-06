using System;
using static System.Math;

public class main{
	static (double,int) integrate
	(Func<double,double> f, double a, double b,
	double δ=0.001, double ε=0.001, double f2=Double.NaN, double f3=Double.NaN)
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
	static double erf(double z){
		if(z<0) return -erf(-z);
		else if(0<=z && z<=1) {
			Func<double,double> integrand = x => Math.Exp(-Math.Pow(x,2));
			(double res, int c) = integrate(integrand,0,z);
			return (2.0/Math.Sqrt(Math.PI))*res;
		}
		else{
			Func<double,double> integrand = t => Math.Exp(-Math.Pow(z+(1-t)/t,2))/t/t;
			(double res, int c) = integrate(integrand,0,1);
			return 1-(2.0/Math.Sqrt(Math.PI))*res;
		}
	}
	static (double,int) ccIntegrate
	(Func<double,double> f, double a, double b,
	double δ=0.001, double ε=0.001)
	{
		Func<double,double> newfunc = x => f( (a+b)/2+(b-a)/2*Math.Cos(x) )*Math.Sin(x)*(b-a)/2;
		return integrate(newfunc,0,Math.PI);
	}

	public static int Main(){
		System.Console.WriteLine("Part A:\n");
		(double sqrtXIntegral, int counts1) = integrate(Math.Sqrt,0,1);
		System.Console.Write("Integral of Sqrt(x) from 0 to 1: ");
		System.Console.Write(sqrtXIntegral);
		System.Console.WriteLine($" which is approximately 2/3 ({counts1} evaluations)");

		System.Console.Write("\n");
		
		Func<double,double> sqrtX1Over = x => 1.0/Math.Sqrt(x);
		(double sqrtX1OverIntegral, int counts2) = integrate(sqrtX1Over,0,1);
		System.Console.Write("Integral of 1/Sqrt(x) from 0 to 1: ");
		System.Console.Write(sqrtX1OverIntegral);
		System.Console.WriteLine($" which is approximately 2 ({counts2} evaluations)");

		System.Console.Write("\n");
		
		Func<double,double> equation3 = x => 4*Math.Sqrt((1-Math.Pow(x,2)));
		(double equation3Integral, int counts3) = integrate(equation3,0,1);
		System.Console.Write("Integral of 4*sqrt(1-x^2) from 0 to 1: ");
		System.Console.Write(equation3Integral);
		System.Console.WriteLine($" which is approximately pi ({counts3} evaluations)");

		System.Console.Write("\n");
		
		Func<double,double> equation4 = x => Math.Log(x)/Math.Sqrt(x);
		(double equation4Integral, int counts4) = integrate(equation4,0,1);
		System.Console.Write("Integral of Ln(x)/Sqrt(x) from 0 to 1: ");
		System.Console.Write(equation4Integral);
		System.Console.WriteLine($" which is approximately -4 ({counts4} evaluations)");

		//We output the values of the found error functions to the error output
		for(double x=-3;x<=3;x+=1.0/8){//we use doubles that can be represented by binary numbers
			System.Console.Error.WriteLine($"{x} {erf(x)}");
		}

		System.Console.WriteLine("\n\nPart B:\n");
		System.Console.WriteLine("We now try the Clenshaw-Curtis variable transformation and compare the number of evaluations of the integrand with the scipy quad function in python:");
		
		(double sqrtX1OverIntegralCC, int counts1CC) = ccIntegrate(sqrtX1Over,0,1);
		System.Console.Write("Integral of 1/Sqrt(x) from 0 to 1: ");
		System.Console.Write(sqrtX1OverIntegralCC);
		System.Console.WriteLine($" which is approximately 2 ({counts1CC} evaluations - 231 using scipy)");

		System.Console.Write("\n");
		
		(double equation4IntegralCC, int counts2CC) = ccIntegrate(equation4,0,1);
		System.Console.Write("Integral of Ln(x)/Sqrt(x) from 0 to 1: ");
		System.Console.Write(equation4IntegralCC);
		System.Console.WriteLine($" which is approximately -4 ({counts2CC} evaluations - 315 using scipy)");

		System.Console.Write("\n");
		System.Console.WriteLine("We see that the number of evaluations is reduced drastically while remaining almost as accurate as before.");
		return 0;
	}
}
