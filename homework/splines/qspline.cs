using System;
using static System.Math;

public class qspline {
	vector x,y,b,c;
	public qspline(vector xs,vector ys){
		/* x=xs.copy(); y=ys.copy(); calculate b and c */
		x=xs.copy();
		y=ys.copy();

		b=new vector(x.size-1);
		c=new vector(x.size-1);
		
		//Forward-recursion:
		c[0] = 0;
		for(int i=0;i<c.size-1;i++){
			double dxi = (x[i+1]-x[i]);
			double dxip1 = (x[i+2]-x[i+1]);
			
			double dyi = (y[i+1]-y[i]);
			double dyip1 = (y[i+2]-y[i+1]);

			double pi = dyi/dxi;
			double pip1 = dyip1/dxip1;
			
			c[i+1]=(1.0/dxip1)*((pip1-pi-c[i]*dxi));
		}

		//Backward-recursion:
		c[c.size-1] /= 2.0;
		for(int i=c.size-2;i>1;i--){
			double dxi = (x[i+1]-x[i]);
			double dxip1 = (x[i+2]-x[i+1]);
			
			double dyi = (y[i+1]-y[i]);
			double dyip1 = (y[i+2]-y[i+1]);

			double pi = dyi/dxi;
			double pip1 = dyip1/dxip1;
			
			c[i] = (1.0/dxi)*(pip1-pi-c[i+1]*dxip1);
		}
		
		for(int i=0;i<b.size;i++) {
			double dxi = (x[i+1]-x[i]);
			double dyi = (y[i+1]-y[i]);
			double pi = dyi/dxi;
			
			b[i] = pi-c[i]*dxi;
		}	
	}//init

	public static int binsearch(vector x, double z)
	{/* locates the interval for z by bisection */
		if( z<x[0] || z>x[x.size-1] ) throw new Exception("binsearch: bad z");
		int i=0, j=x.size-1;
		while(j-i>1){
			int mid=(i+j)/2;
			if(z>x[mid]) i=mid; else j=mid;
			}
		return i;
	}//binsearch

	public double evaluate(double z){
		/* evaluate the spline */
		
		int i = binsearch(x,z);	
		
		double si = y[i]+b[i]*(z-x[i])+c[i]*Pow(z-x[i],2.0);
		
		return si;
	}//evaluate
	
	public double derivative(double z){
		/* evaluate the derivative */
		//Derivative of si(x) is bi+2*ci*(x-xi)

		int i = binsearch(x,z);
		double dsi = b[i]+2*c[i]*(z-x[i]);
		return dsi;
	}//derivative

	public double integral(double z){
		/* evaluate the integral */
		int i = binsearch(x,z);

		double intOut = 0;
		for(int j=0;j<i;j++){  
		    intOut += (y[j]*(x[j+1]-x[j])+b[j]*Pow((x[j+1]-x[j]),2)/2+c[j]*Pow((x[j+1]-x[j]),3)/3);
		}

		//Need to include point at z:
		intOut += (y[i]*(z-x[i])+b[i]*Pow((z-x[i]),2)/2+c[i]*Pow((z-x[i]),3)/3);

		//integration gives us an additional constant which in this case is C = -1.0:
		double C = -1.0;
		intOut += C;

		return intOut;
		
	}//integral


public static int Main(){
		Random rnd = new Random(123);
		vector xs = new vector(11);	
		vector ys = new vector(11);
	
		Func<double,double> test_function = x => Sin(x);
		
		for(int i=0;i<xs.size;i++){
			xs[i]=i;
			ys[i]=test_function(xs[i]);
			Console.WriteLine($"{xs[i]} {ys[i]}");
		}
		
		System.Console.Out.WriteLine("\n");

		int numPoints = 100;
		
		double start = xs[0];
		double last = xs[xs.size-1];
		double step = (last-start)/(numPoints-1);
		
		//Sin output
		for(int i=0;i<numPoints;i++){
			double usex = start + i*step;
			System.Console.Out.WriteLine($"{usex} {test_function(usex)}");
		}

		System.Console.Out.WriteLine("\n");

		//Quadratic interp output
		var qspline = new qspline(xs,ys);
		
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			System.Console.Out.WriteLine($"{usex} {qspline.evaluate(usex)}");
		}

		System.Console.Out.WriteLine("\n");
	
		//Derivative output
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			System.Console.Out.WriteLine($"{usex} {qspline.derivative(usex)}");
		}
		
		System.Console.Out.WriteLine("\n");
	
		//Integral output
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			System.Console.Out.WriteLine($"{usex} {qspline.integral(usex)}");
		}


		return 0;
	}
}

