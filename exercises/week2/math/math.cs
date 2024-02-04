using System;
using static System.Console;
using static System.Math;

class math{
	static void Main(){
		double sqrt2=Sqrt(2.0);
		Write($"sqrt2^2 = {sqrt2*sqrt2} (should equal 2)\n");

		double power215 = Pow(2.0,1.0/5.0);
		Write($"power of 2 to 1/5 = {power215:F5} (should equal 1.148...)\n");

		double epi = Exp(PI);
		Write($"power of e to pi = {epi:F5} should be slightly larger than e*e*e = {Exp(3):F5}\n");
		
		double pie = Pow(PI,E);
		Write($"power of pi to e = {pie:F5} should be slightly smaller than pi*pi*pi = {(PI*PI*PI):F5}\n");

		for(int i=1;i<=10;i++){
			Write($"For {i}:\n");
			Write($"G({i}) = {sfuns.fgamma(i):F6}\n");
			Write($"lnG({i}) = {sfuns.flngamma(i):F1}, giving exp(lnG({i})) = {Exp(sfuns.flngamma(i)):F1}\n");
		}
	}
}
