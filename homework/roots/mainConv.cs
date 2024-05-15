using System;
using static System.Math;

public class mainConv{
	public static int Main(string[] args) {
		double rmin = 0;
		double rmax = 0;
		double acc = 0;
		double eps = 0;

		foreach(var arg in args){
			var words = arg.Split(':');
			if(words[0]=="-rmin"){
				rmin = double.Parse(words[1]);
			}
			else if(words[0]=="-rmax"){
				rmax = double.Parse(words[1]);
			}
			else if(words[0]=="-acc"){
				acc = double.Parse(words[1]);
			}
			else if(words[0]=="-eps"){
				eps = double.Parse(words[1]);
			}
		}
		
		vector odestart = new vector(rmin-rmin*rmin,1-2*rmin);

		Func<vector, vector> M_E = delegate(vector E) {
			Func<double, vector, vector> f_E = delegate(double x, vector y){
				//Just like in the ODE task we can rewrite the second order system to an ODE by
				//introducing y0=f and y1=f'
				//Our output will now become (f,f')
				return new vector(y[1],-2*(E[0]+(1.0/x))*y[0]);
			};
			
			var (_xlist,_ylist) = ODE.driver(f_E,(rmin,rmax),odestart,acc=acc,eps=eps);
			
			vector res = new vector(1);
			res[0] = _ylist[_ylist.size-1][0];
			return res;
		};
		
		vector startE = new vector("-5");

		vector resultE = roots.newton(M_E,startE);

		Console.WriteLine($"{rmin} {rmax} {acc} {eps} {resultE[0]}");

		return 0;
	}
}
