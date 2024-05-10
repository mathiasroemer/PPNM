using System;
using static System.Math;

public class linspline {
	public static double linterp(double[] x, double[] y, double z){
		int i=binsearch(x,z);
		double dx=x[i+1]-x[i]; if(!(dx>0)) throw new Exception("uups...");
		double dy=y[i+1]-y[i];
		return y[i]+dy/dx*(z-x[i]);
	}

	public static int binsearch(double[] x, double z)
	{/* locates the interval for z by bisection */
		if( z<x[0] || z>x[x.Length-1] ) throw new Exception("binsearch: bad z");
		int i=0, j=x.Length-1;
		while(j-i>1){
			int mid=(i+j)/2;
			if(z>x[mid]) i=mid; else j=mid;
			}
		return i;
	}

	public static double linterpInteg(double[] x, double[] y, double z){
		//integral of y=ax+b from x0 to z is -(1/2)*(x_0-z)*(a*(x_0+z)+2*b)
		
		int i = binsearch(x,z);	
		double intOut = 0;
		
		//Sum of all the integral parts:
		for(int j=0;j<i;j++){
			double a = (y[j+1]-y[j])/(x[j+1]-x[j]);
			double b = y[j]-a*x[j];
			intOut -= (1.0/2)*(x[j]-x[j+1])*(a*(x[j]+x[j+1])+2*b);
		}

		//We need to include part at z:
		double anew = (y[i+1]-y[i])/(x[i+1]-x[i]);
		double bnew = y[i]-anew*x[i];
		intOut -= (1.0/2)*(x[i]-z)*(anew*(x[i]+z)+2*bnew);


		return intOut;
	}
	
	public static int Main(){
		double[] x = {0,1,2,3,4,5,6,7,8,9};
		double[] y = new double[x.Length];
		for(int i=0;i<x.Length;i++)y[i]=Cos(x[i]);
	
		//Data output
		for(int i=0;i<x.Length;i++){
			System.Console.Out.WriteLine($"{x[i]} {y[i]}");
		}

		System.Console.Out.WriteLine("\n");
		
		int numPoints = 100;
		
		double start = x[0];
		double last = x[x.Length-1];
		double step = (last-start)/(numPoints-1);
		
		//Cos output
		for(int i=0;i<numPoints;i++){
			double usex = start + i*step;
			System.Console.Out.WriteLine($"{usex} {Cos(usex)}");
		}

		System.Console.Out.WriteLine("\n");

		//Linterp output
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			System.Console.Out.WriteLine($"{usex} {linterp(x,y,usex)}");
		}

		System.Console.Out.WriteLine("\n");

		//Anti-derivative output
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			System.Console.Out.WriteLine($"{usex} {linterpInteg(x,y,usex)}");
		}

		return 0;
	}
}
