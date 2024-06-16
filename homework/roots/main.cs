using System;
using static System.Math;

public class main{
	public static int Main(){
		Console.WriteLine("Part A:\n");
		Console.WriteLine("The Rosenbrock's valley function is: f(x,y) = (1-x)^2 + 100*(y-x^2)^2");
		Console.WriteLine("The gradient is the found analytically to be ∇f(x,y) = (2*(-1+x+200*x^3 - 200*x*y),200*(-x^2+y))");
		Func<vector,vector> GradRbvFunc = delegate(vector x) {
			var _x = x[0];
			var _y = x[1];
			return new vector(2*(-1+_x+200*Math.Pow(_x,3)-200*_x*_y),200*(-Math.Pow(_x,2)+_y));
		};

		vector start = new vector(-10,10);
	
		vector newt = roots.newton(GradRbvFunc, start);

		newt.print("Using Newton's method with a starting guess at (-10,10) we find the global minimum at:\n");
		Console.WriteLine("in agreement with the provided minimum at (a,a^2) given by the wikipedia.");
		Console.WriteLine("\n");
		Console.WriteLine("We now want to use Newton's method to find the minima of the Himmelblau's function: f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2");
		Console.WriteLine("Which has an analytic gradient of: ∇f(x,y) = (2*(2*x*(x^2+y-11)+x+y^2-7),2*(x^2+2*y*(x+y^2-7)+y-11))");

		Func<vector,vector> GradHimmelFunc = delegate(vector x) {
			var _x = x[0];
			var _y = x[1];
			return new vector(2*(2*_x*(Math.Pow(_x,2)+_y-11)+_x+Math.Pow(_y,2)-7),2*(Math.Pow(_x,2)+2*_y*(_x+Math.Pow(_y,2)-7)+_y-11));
		};

		int[] startX = {10,-7,-12,9};
		int[] startY = {10,7,-12,-9};
		
		Console.WriteLine("We try now four different intial guesses and get the results:");

		for(int i = 0; i<startX.Length; i++){
			Console.Write("\n");
			Console.WriteLine($"Initial guess: ({startX[i]}, {startY[i]})");
			Console.WriteLine("Result:");

			vector startGuess = new vector(startX[i],startY[i]);
			
			var result = roots.newton(GradHimmelFunc,startGuess);

			result.print("\n");
		}

		Console.Write("\n");
		Console.WriteLine("Which are all in agreement with the minima given by the wikipedia.\n");

		Console.WriteLine("Part B:\n");
		Console.WriteLine("We now want to find the lowest root of M(E) = F_E(rmax) where we find F_E by using our ODE solver on the s-wave radial Schrödinger equation for hydrogen: -(1/2)f'' - (1/r)f = Ef\n");
		double rmin = 0.1;
		double rmax = 8.0;
		vector odestart = new vector(rmin-rmin*rmin,1-2*rmin);

		Func<vector, vector> M_E = delegate(vector E) {
			Func<double, vector, vector> f_E = delegate(double x, vector y){
				//Just like in the ODE task we can rewrite the second order system to an ODE by
				//introducing y0=f and y1=f'
				//Our output will now become (f,f')
				return new vector(y[1],-2*(E[0]+(1.0/x))*y[0]);
			};
			
			var (_xlist,_ylist) = ODE.driver(f_E,(rmin,rmax),odestart);
			
			vector res = new vector(1);
			res[0] = _ylist[_ylist.size-1][0];
			return res;
		};
		
		vector startE = new vector("-5");

		Console.WriteLine($"We first we guess that the root energy is at {startE[0]} (a better guess might be -1/2 but we don't know that yet)");

		vector resultE = roots.newton(M_E,startE);
		
		Console.WriteLine($"But then after using our routine find that the lowest root energy E_0 is {resultE[0]:f3}");

		Func<double, vector, vector> F_E = delegate(double x, vector y){
			return new vector(y[1],-2*(resultE[0]+(1.0/x))*y[0]);
		};
			
		var (xlist,ylist) = ODE.driver(F_E,(rmin,rmax),odestart);
	
		Console.WriteLine("We then plot the result wavefunction for this found root energy together with the exact result wave function");

		for(int i=0;i<xlist.size;i++){
			Console.Error.WriteLine($"{xlist[i]} {ylist[i][0]}");
		}

		Console.WriteLine("\nPart C:");

		Console.WriteLine("We now try to use the Newton's method with quadratic line-search to find the minima of the Himmelblau's function: f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2");
		for(int i = 0; i<startX.Length; i++){
			Console.Write("\n");
			Console.WriteLine($"Initial guess: ({startX[i]}, {startY[i]})");
			Console.WriteLine("Result:");

			vector startGuess = new vector(startX[i],startY[i]);
			
			var result = roots.newtonQuad(GradHimmelFunc,startGuess);

			result.print("\n");
		}
	
		Console.Write("\n");
		Console.WriteLine("Which are all in agreement with what we found earlier.\n");

	return 0;
	}
}
