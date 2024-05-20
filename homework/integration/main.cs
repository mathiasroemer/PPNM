using System;
using static System.Math;

public class main{
	static (double,double,int) integrate
	(Func<double,double> f, double a, double b,
	double δ=0.001, double ε=0.001, double f2=Double.NaN, double f3=Double.NaN)
	{
		int evals = 0;

		if(a == double.NegativeInfinity && b == double.PositiveInfinity){
			Func<double,double> newfunc = t => f(t/(1-Pow(t,2)))*((1+Pow(t,2))/(Pow((1-Pow(t,2)),2)));
			return integrate(newfunc,-1,1);
		}
		else if(a == double.NegativeInfinity) {
			Func<double,double> newfunc = t => f(b-((1-t)/t))*(1.0/Pow(t,2));
			return integrate(newfunc,0,1);		
		}
		else if(b == double.PositiveInfinity) {
			Func<double,double> newfunc = t => f((a+(1-t)/t))*(1.0/(Pow(t,2)));
			return integrate(newfunc,0,1);
		}

		double h=b-a;
		if(Double.IsNaN(f2)){ f2=f(a+2*h/6); f3=f(a+4*h/6); evals+=2;} // first call, no points to reuse
		double f1=f(a+h/6), f4=f(a+5*h/6); evals+=2;
		double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // higher order rule
		double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule
		double err = Math.Abs(Q-q);
		if (err <= δ+ε*Math.Abs(Q)) return (Q,err,evals);
		else {
			(double int1, double error1, int evals1) = integrate(f,a,(a+b)/2,δ/Math.Sqrt(2),ε,f1,f2);
			(double int2, double error2, int evals2) = integrate(f,(a+b)/2,b,δ/Math.Sqrt(2),ε,f3,f4);
			return (int1+int2,Sqrt(Pow(error1,2)+Pow(error2,2)),evals+evals1+evals2);
		}
	}
	static double erf(double z){
		if(z<0) return -erf(-z);
		else if(0<=z && z<=1) {
			Func<double,double> integrand = x => Math.Exp(-Math.Pow(x,2));
			(double res, double error, int c) = integrate(integrand,0,z);
			return (2.0/Math.Sqrt(Math.PI))*res;
		}
		else{
			Func<double,double> integrand = t => Math.Exp(-Math.Pow(z+(1-t)/t,2))/t/t;
			(double res, double error, int c) = integrate(integrand,0,1);
			return 1-(2.0/Math.Sqrt(Math.PI))*res;
		}
	}
	static (double,double,int) ccIntegrate
	(Func<double,double> f, double a, double b,
	double δ=0.001, double ε=0.001)
	{
		Func<double,double> newfunc = x => f( (a+b)/2+(b-a)/2*Math.Cos(x) )*Math.Sin(x)*(b-a)/2;
		return integrate(newfunc,0,Math.PI);
	}

	public static int Main(){
		System.Console.WriteLine("Part A:\n");
		(double sqrtXIntegral, double sqrtXIntegral_error, int counts1) = integrate(Math.Sqrt,0,1);
		System.Console.Write("Integral of Sqrt(x) from 0 to 1: ");
		System.Console.Write($"{sqrtXIntegral:F4} +- {sqrtXIntegral_error:F4}");
		System.Console.WriteLine($" which is approximately 2/3 ({counts1} evaluations)");

		System.Console.Write("\n");
		
		Func<double,double> sqrtX1Over = x => 1.0/Math.Sqrt(x);
		(double sqrtX1OverIntegral, double sqrtX1OverIntegral_error, int counts2) = integrate(sqrtX1Over,0,1);
		System.Console.Write("Integral of 1/Sqrt(x) from 0 to 1: ");
		System.Console.Write($"{sqrtX1OverIntegral:F4} +- {sqrtX1OverIntegral_error:F4}");
		System.Console.WriteLine($" which is approximately 2 ({counts2} evaluations)");

		System.Console.Write("\n");
		
		Func<double,double> equation3 = x => 4*Math.Sqrt((1-Math.Pow(x,2)));
		(double equation3Integral, double equation3Integral_error, int counts3) = integrate(equation3,0,1);
		System.Console.Write("Integral of 4*sqrt(1-x^2) from 0 to 1: ");
		System.Console.Write($"{equation3Integral:F4} +- {equation3Integral_error:F4}");
		System.Console.WriteLine($" which is approximately pi ({counts3} evaluations)");

		System.Console.Write("\n");
		
		Func<double,double> equation4 = x => Math.Log(x)/Math.Sqrt(x);
		(double equation4Integral, double equation4Integral_error, int counts4) = integrate(equation4,0,1);
		System.Console.Write("Integral of Ln(x)/Sqrt(x) from 0 to 1: ");
		System.Console.Write($"{equation4Integral:F4} +- {equation4Integral_error:F4}");
		System.Console.WriteLine($" which is approximately -4 ({counts4} evaluations)");

		//We output the values of the found error functions to the error output
		for(double x=-3;x<=3;x+=1.0/8){//we use doubles that can be represented by binary numbers
			System.Console.Error.WriteLine($"{x} {erf(x)}");
		}

		System.Console.WriteLine("\n\nPart B:\n");
		System.Console.WriteLine("We now try the Clenshaw-Curtis variable transformation and compare the number of evaluations of the integrand with the scipy quad function in python:");
		
		(double sqrtX1OverIntegralCC, double sqrtX1OverIntegralCC_error, int counts1CC) = ccIntegrate(sqrtX1Over,0,1);
		System.Console.Write("Integral of 1/Sqrt(x) from 0 to 1: ");
		System.Console.Write($"{sqrtX1OverIntegralCC:F4} +- {sqrtX1OverIntegralCC_error:F4}");
		System.Console.WriteLine($" which is approximately 2 ({counts1CC} evaluations - 231 using scipy)");

		System.Console.Write("\n");
		
		(double equation4IntegralCC, double equation4IntegralCC_error, int counts2CC) = ccIntegrate(equation4,0,1);
		System.Console.Write("Integral of Ln(x)/Sqrt(x) from 0 to 1: ");
		System.Console.Write($"{equation4IntegralCC:F4} +- {equation4IntegralCC_error:F4}");
		System.Console.WriteLine($" which is approximately -4 ({counts2CC} evaluations - 315 using scipy)");

		System.Console.Write("\n");
		System.Console.WriteLine("We see that the number of evaluations is reduced drastically while remaining almost as accurate as before.");
		System.Console.WriteLine("\n\nPart C:\n");
		System.Console.WriteLine("For the error part, the integration errors has been distributed out into Part A and Part B\n");
		System.Console.WriteLine("Now for the infinite limit integrals:\n");
		System.Console.WriteLine("We first test the case where both limits a infinity on the gaussian integral:");

		Func<double,double> gaussian = x => Exp(-x*x);
		(double gaussianIntegral, double gaussianIntegral_error, int counts5) = integrate(gaussian,double.NegativeInfinity,double.PositiveInfinity);
		System.Console.Write("Integral of Exp(-x^2) from -Inf to +Inf: ");
		System.Console.Write($"{gaussianIntegral:F4} +- {gaussianIntegral_error:F4}");
		System.Console.WriteLine($" which is approximately Sqrt(PI) = {Sqrt(PI):F4} which is the expected result ({counts5} evaluations)");
		System.Console.WriteLine("Using scipy integration we get a result of 1.7725 with 270 evaluations.\n");

		System.Console.WriteLine("We now test the case where only the upper limit is infinity on the exponential decay integral:");
		System.Console.Write("Integral of Exp(-x) from 0 to +Inf: ");
	
		Func<double,double> expDecay = x => Exp(-x);
		(double expDecayIntegral, double expDecayIntegral_error, int counts6) = integrate(expDecay,0,double.PositiveInfinity);
		System.Console.Write($"{expDecayIntegral:F4} +- {expDecayIntegral_error:F4}");
		System.Console.WriteLine($" which is approximately 1 which is the expected result ({counts6} evaluations)");
		System.Console.WriteLine("Using scipy integration we get a result of 1 with 135 evaluations.\n");
		
		System.Console.WriteLine("At last we test the case where only the lower limit is infinity on the exponential decay integral with finite upper bound:");
		System.Console.Write("Integral of Exp(x) from -Inf to 0: ");
	
		Func<double,double> expDecayFUB = x => Exp(x);
		(double expDecayFUBIntegral, double expDecayFUBIntegral_error, int counts7) = integrate(expDecayFUB,double.NegativeInfinity,0);
		System.Console.Write($"{expDecayFUBIntegral:F4} +- {expDecayFUBIntegral_error:F4}");
		System.Console.WriteLine($" which is approximately 1 which is again the expected result ({counts7} evaluations)");
		System.Console.WriteLine("Using scipy integration we get a result of 1 with 135 evaluations.\n");

		return 0;
	}
}
