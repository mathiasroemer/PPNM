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

	public static vector gradient_central(Func<vector,double> phi,vector x){
		vector gphi = new vector(x.size);
		vector xp = x.copy();
		vector xm = x.copy();

		for(int i=0;i<x.size;i++){
			double dx=Max(Abs(x[i]),1)*Pow(2,-26);
			xp[i] += dx;
			xm[i] -= dx;

			gphi[i] = (phi(xp)-phi(xm))/(2*dx);

			xp[i] -= dx;
			xm[i] += dx;
		}
		return gphi;
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

	public static matrix hessian_central(Func<vector,double> phi,vector x){
		matrix H=new matrix(x.size);
		vector gphix=gradient_central(phi,x);
		double phix = phi(x);
	
		vector xpp = x.copy();
		vector xmm = x.copy();
		vector xmp = x.copy();
		vector xpm = x.copy();

		for(int j=0;j<x.size;j++){
			for(int i=0;i<x.size;i++) { 
				double dxj=Max(Abs(x[j]),1)*Pow(2,-13); /* for numerical gradient */
				double dxi=Max(Abs(x[i]),1)*Pow(2,-13);
		
				xpp[i] += dxi;
				xpp[j] += dxj;
				
				xpm[i] += dxi;
				xpm[j] -= dxj;

				xmp[i] -= dxi;
				xmp[j] += dxj;

				xmm[i] -= dxi;
				xmm[j] -= dxj;	

				H[i,j] = (phi(xpp)-phi(xpm)-phi(xmp)+phi(xmm))/(4.0*dxi*dxj);

				xpp[i] -= dxi;
				xpp[j] -= dxj;
				
				xpm[i] -= dxi;
				xpm[j] += dxj;

				xmp[i] += dxi;
				xmp[j] -= dxj;

				xmm[i] += dxi;
				xmm[j] += dxj;

			}
		}
		//return H;
		return (H+H.T)/2; // you think?
	}

	public static (vector,int) newton(
		Func<vector,double> φ, /* objective function */
		vector x,              /* starting point */
		bool central=false,
		double acc=1e-3        /* accuracy goal, on exit |∇φ| should be < acc */
	){
		int NSteps = 0;
		int maxsteps = 1000;
		do{ /* Newton's iterations */
			NSteps++;
			vector gφ = new vector(x.size);
			if(central==true){
				gφ = gradient_central(φ,x);
			}
			else {
				gφ = gradient(φ,x);
			}
			if(gφ.norm() < acc) break; /* job done */
			matrix H=new matrix(x.size);
			if(central==true) {
				H = hessian_central(φ,x);	
			}
			else {
				H = hessian(φ,x);
			
			}
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
		
		int maxSteps = 0;
		for(int i=0; i<startXRbv.Length;i++) {
			vector start = new vector(startXRbv[i],startYRbv[i]);
	
			(vector newt, int steps) = newton(RbvFunc, start);
			
			newt.print($"Start guess: ({startXRbv[i]}, {startYRbv[i]}), Global minimum found after {steps} steps:\n");
			if(steps>maxSteps) maxSteps = steps;
		}
		Console.WriteLine($"\nWe see that the minimum is closely around the provided minimum at (a,a^2) = (1,1) given by the wikipedia and that the algorithm stops after {maxSteps} steps in case the convergence criterion cannot be reached.");
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


		Console.WriteLine("\nPart B:\n");
		Console.WriteLine("We want to fit the Higgs Boson data to the Breit-Wigner function\n");
		var energies = new genlist<double>();
		var signals = new genlist<double>();
		var errors  = new genlist<double>();
		var separators = new char[] {' ','\t'};
		var options = StringSplitOptions.RemoveEmptyEntries;
		do{
			string line=Console.In.ReadLine();
			if(line==null)break;
			string[] words=line.Split(separators,options);
			energies.add(double.Parse(words[0]));
			signals.add(double.Parse(words[1]));
			errors.add(double.Parse(words[2]));
		}while(true);
	
		Func<double,vector,double> Breit_Wigner = delegate(double energy, vector param) {
			double E = energy;
			double m = param[0];
			double G = param[1];
			double A = param[2];
			return A/(Pow(E-m,2)+Pow(G,2.0)/4.0);
		};

		Func<vector,double> deviation_func = delegate(vector param) {
			double sum = 0.0;
			for(int i=0;i<energies.size;i++){
				double Ei = energies[i];
				double si = signals[i];
				double dsi = errors[i];

				sum += Pow((Breit_Wigner(Ei,param)-si)/dsi,2.0);
			}
			return sum;
		};
		
		vector BW_start = new vector(126,2,5.5); //mass, width, scale-factor

		(vector BW_result, int BW_steps) = newton(deviation_func,BW_start,false,0.00001);
	
		BW_start.print("Start parameters:");
		Console.WriteLine("(Should be a fair guess looking at the data)\n");
		Console.WriteLine($"After {BW_steps} minimization steps we find the final parameters:");
		BW_result.print();
		
		Console.WriteLine($"\nWe find a mass of {BW_result[0]:F2} GeV with a resonance width of {BW_result[1]:F2} and a scale factor of {BW_result[2]:F2}");

		for(double i=91;i<=170;i+=1.0/8){
			Console.Error.WriteLine($"{i} {Breit_Wigner(i,BW_result)}");
		}
		
		Console.WriteLine("\nPart C:\n");
		Console.WriteLine("We now do the same as for Part A, just with a central finite difference instead of a forward difference");

		Console.WriteLine("The Rosenbrock's valley function is: f(x,y) = (1-x)^2 + 100*(y-x^2)^2\n");

		for(int i=0; i<startXRbv.Length;i++) {
			vector start = new vector(startXRbv[i],startYRbv[i]);
	
			(vector newt, int steps) = newton(RbvFunc, start,true);
			
			newt.print($"Start guess: ({startXRbv[i]}, {startYRbv[i]}), Global minimum found after {steps} steps:\n");
			if(steps>maxSteps) maxSteps = steps;
		}
		Console.WriteLine($"\nWe see that the minimum is closely around the provided minimum at (a,a^2) = (1,1) given by the wikipedia and that it this time got the minimum of all the start guesses.");
		Console.WriteLine("\n");
		Console.WriteLine("We now want to use Newton's method with a central gradient to find the minima of the Himmelblau's function: f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2");
		
		Console.WriteLine("We try now four different intial guesses and get the results:");

		for(int i = 0; i<startXHB.Length; i++){
			Console.Write("\n");
			Console.WriteLine($"Initial guess: ({startXHB[i]}, {startYHB[i]})");

			vector startGuess = new vector(startXHB[i],startYHB[i]);
			
			(vector result, int steps) = newton(HimmelFunc,startGuess,true);
			
			Console.WriteLine($"Result after {steps} steps:");

			result.print("\n");
		}

		Console.Write("\n");
		Console.WriteLine("Which are all in agreement with the minima we go from the forward difference.\n");
		Console.WriteLine("We see that the central finite difference approximations for the derivatives is better than the forward difference approximations since we not get results for all our guesses.");

		return 0;
	}
}
