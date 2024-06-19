using System;
using static System.Math;

public class main {

public static int Main(){
		vector xs = new vector(-3,-2,-1,-0.25,0.25,1,2,3);

		vector ys = new vector(xs.size);
		vector dys = new vector(xs.size);

		Func<double,double> sign_function = x => Math.Sign(x);
		Func<double,double> sign_function_derivative = x => 0;

		for(int i=0;i<xs.size;i++){
			ys[i]=sign_function(xs[i]);
			dys[i]=sign_function_derivative(xs[i]);
			Console.WriteLine($"{xs[i]} {ys[i]}");
		}
		
		System.Console.Out.WriteLine("\n");

		int numPoints = 100;
		
		double start = xs[0];
		double last = xs[xs.size-1];
		double step = (last-start)/(numPoints-1);
		
		//Sign output
		for(int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {sign_function(usex)}");
		}

		System.Console.Out.WriteLine("\n");
		///Subspline interp output on sign
		var subspline = new subspline(xs,ys,dys);
		
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {subspline.evaluate(usex)}");
		}
		
		System.Console.Out.WriteLine("\n");
		///cspline interp output on sign
		var cspline = new cspline(xs,ys);
		
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {cspline.evaluate(usex)}");
		}

		System.Console.Out.WriteLine("\n");
		
		//Test the derivative and integral on the Cos function
		int NPoints = 15;
		vector xsCos = new vector(NPoints);

		vector ysCos = new vector(xsCos.size);
		vector dysCos = new vector(xsCos.size);

		Func<double,double> cos_function = x => Math.Cos(x);
		Func<double,double> cos_function_derivative = x => -Math.Sin(x);

		for(int i=0;i<xsCos.size;i++){
			xsCos[i]=i;
			ysCos[i]=cos_function(xsCos[i]);
			dysCos[i]=cos_function_derivative(xsCos[i]);
			Console.WriteLine($"{xsCos[i]} {ysCos[i]}");
		}
		
		System.Console.Out.WriteLine("\n");

		
		start = xsCos[0];
		last = xsCos[xsCos.size-1];
		step = (last-start)/(numPoints-1);
		
		//Cos output
		for(int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {cos_function(usex)}");
		}

		System.Console.Out.WriteLine("\n");

		///Subspline interp output on cos
		var subspline_cos = new subspline(xsCos,ysCos,dysCos);
			
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {subspline_cos.evaluate(usex)}");
		}
		
		System.Console.Out.WriteLine("\n");
			
		//Subspline derivative
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {subspline_cos.derivative(usex)}");
		}
		
		System.Console.Out.WriteLine("\n");
		
		//Subspline Integral
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {subspline_cos.integral(usex)}");
		}
		
		System.Console.Out.WriteLine("\n");
		
		return 0;
	}
}

