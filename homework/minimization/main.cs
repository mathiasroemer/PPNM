using System;
using static System.Math;

public class main{
	public static vector gradient(Func<vector,double> φ,vector x){
		vector gφ = new vector(x.size);
		double φx = φ(x); /* no need to recalculate at each step */
		for(int i=0;i<x.size;i++){
			double dx=Max(Abs(x[i]),1)*Pow(2,-26);
			x[i]+=dx;
			gφ[i]=(φ(x)-φx)/dx;
			x[i]-=dx;
		}
		return gφ;
	}

	public static matrix hessian(Func<vector,double> φ,vector x){
		matrix H=new matrix(x.size);
		vector gφx=gradient(φ,x);
		for(int j=0;j<x.size;j++){
			double dx=Max(Abs(x[j]),1)*Pow(2,-13); /* for numerical gradient */
			x[j]+=dx;
			vector dgφ=gradient(φ,x)-gφx;
			for(int i=0;i<x.size;i++) H[i,j]=dgφ[i]/dx;
			x[j]-=dx;
		}
		//return H;
		return (H+H.T)/2; // you think?
	}

	public static (vector,int) newton(
		Func<vector,double> φ, /* objective function */
		vector x,              /* starting point */
		double acc=1e-3        /* accuracy goal, on exit |∇φ| should be < acc */
	){
		int NSteps = 0;
		int maxsteps = 1000;
		do{ /* Newton's iterations */
			NSteps++;
			var gφ = gradient(φ,x);
			if(gφ.norm() < acc) break; /* job done */
			var H = hessian(φ,x);
			//var QRH = givensQR(H);   /* QR decomposition */
			//var dx = QRH.solve(-∇φ); /* Newton's step */
			var QRH = QRGS.decomp(H);
			var (Q,R) = QRH;
			vector dx = QRGS.solve(Q,R,-gφ);
		
			double λ=1,φx=φ(x);
			double λmin = 1.0/1024;
			do{ /* linesearch */
				if( φ(x+λ*dx) < φx ) break; /* good step: accept */
				if( λ < λmin ) break; /* accept anyway */
				λ/=2;
			}while(true);
			x+=λ*dx;
		}while(NSteps < maxsteps);
		return (x,NSteps);
	}//newton
	public static int Main(){
		Console.WriteLine("Part A:\n");
		Console.WriteLine("The Rosenbrock's valley function is: f(x,y) = (1-x)^2 + 100*(y-x^2)^2\n");
		
		Func<vector,double> RbvFunc = delegate(vector x) {
			var _x = x[0];
			var _y = x[1];
			return Math.Pow((1-_x),2) + 100*Math.Pow((_y-Math.Pow(_x,2)),2);
		};

		Console.WriteLine("We now use Newton's method with a numerical gradient in order to find the minimum. This is done using five different initial start guesses:\n"); 
		
		int[] startXRbv = {2,3,-3,4,-10};
		int[] startYRbv = {2,3,3,-3,10};
		
		for(int i=0; i<startXRbv.Length;i++) {
			vector start = new vector(startXRbv[i],startYRbv[i]);
	
			(vector newt, int steps) = newton(RbvFunc, start);
			
			newt.print($"Start guess: ({startXRbv[i]}, {startYRbv[i]}), Global minimum found after {steps} steps:\n");
		}
		Console.WriteLine("\nWe see that the minimum is closely around the provided minimum at (a,a^2) given by the wikipedia and that the algorithm stops after 1000 steps in case the convergence criterion cannot be reached.");
		Console.WriteLine("\n");
		Console.WriteLine("We now want to use Newton's method with a numerical gradient to find the minima of the Himmelblau's function: f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2");
		
		Func<vector,double> HimmelFunc = delegate(vector x) {
			var _x = x[0];
			var _y = x[1];
			return Math.Pow(Math.Pow(_x,2)+_y-11,2)+Math.Pow((_x+Math.Pow(_y,2)-7),2);
		};
		int[] startXHB = {10,-7,-12,9};
		int[] startYHB = {10,7,-12,-9};
		
		Console.WriteLine("We try now four different intial guesses and get the results:");

		for(int i = 0; i<startXHB.Length; i++){
			Console.Write("\n");
			Console.WriteLine($"Initial guess: ({startXHB[i]}, {startYHB[i]})");

			vector startGuess = new vector(startXHB[i],startYHB[i]);
			
			(vector result, int steps) = newton(HimmelFunc,startGuess);
			
			Console.WriteLine($"Result after {steps} steps:");

			result.print("\n");
		}

		Console.Write("\n");
		Console.WriteLine("Which are all in agreement with the minima given by the wikipedia.");


		return 0;
	}
}
