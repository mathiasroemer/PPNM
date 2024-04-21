using System;
using static System.Math;

public class main{
	
	public static matrix jacobian
	(Func<vector,vector> f,vector x,vector fx=null,vector dx=null){
		if(dx == null) dx = x.map(xi => Abs(xi)*Pow(2,-26));
		if(fx == null)fx=f(x);
		matrix J=new matrix(x.size);
		for(int j=0;j < x.size;j++){
			x[j]+=dx[j];
			vector df=f(x)-fx;
			for(int i=0;i < x.size;i++) J[i,j]=df[i]/dx[j];
			x[j]-=dx[j];
			}
		return J;
	}//jacobian
	
	public static vector newton(
	Func<vector,vector>f /* the function to find the root of */
	,vector x            /* the start point */
	,double acc=1e-2     /* accuracy goal: on exit ‖f(x)‖ should be <acc */
	,vector δx=null      /* optional δx-vector for calculation of jacobian */
	){
		vector fx=f(x),z,fz;
		do{ /* Newton's iterations */
			if(fx.norm() < acc) break; /* job done */
			matrix J=jacobian(f,x,fx,δx);
			var QRofJ = QRGS.decomp(J);
			var (Q,R) = QRofJ;
			vector Dx = QRGS.solve(Q,R,-fx);
			double λ=1;
			double λmin=1.0/64; 
			do{ /* linesearch */
				z=x+λ*Dx;
				fz=f(z);
				if( fz.norm() < (1-λ/2)*fx.norm() ) break;
				if( λ < λmin ) break;
				λ/=2;
				}while(true);
			x=z; fx=fz;
			}while(true);
		return x;
	}//newton
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
	
		vector newt = newton(GradRbvFunc, start);

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
			
			var result = newton(GradHimmelFunc,startGuess);

			result.print("\n");
		}

		Console.Write("\n");
		Console.WriteLine("Which are all in agreement with the minima given by the wikipedia.");

		return 0;
	}
}
