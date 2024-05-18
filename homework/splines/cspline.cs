using System;
using static System.Math;

public class cspline {
	vector x,y,b,c,d;
	int n;
	public cspline(vector xs,vector ys){
		n = xs.size;

		x = xs.copy();
		y = ys.copy();
		
		b = new vector(n);
		c = new vector(n-1);
		d = new vector(n-1);

		for ( int i =0;i<n ; i++){x [ i ]=xs[ i ] ; y[ i ]=ys[ i ] ; }

		double[] h = new double[n-1];
		double[] p = new double[n-1]; // VLA
		
		for ( int i =0;i<n-1; i++){h[ i ]=x[i+1]-x[ i ];}
		
		for ( int i =0;i<n-1; i++) p[ i ]=(y[ i+1]-y[ i ] ) / h[ i ] ;
		
		double[] D = new double[n];
		double[] Q = new double[n-1];
	  	double[] B = new double[n]; // building the tridiagonal system :
		
		D[0]=2; 
		for ( int i =0;i<n-2; i++)D[ i +1]=2*h [ i ]/h [ i +1]+2; D[ n-1]=2;
		
		Q[0]=1; 
		
		for ( int i =0;i<n-2; i++) Q[ i+1]=h [ i ]/h [ i +1];
		
		for ( int i =0;i<n-2; i++) B[ i +1]=3*(p[ i ]+p [ i +1]*h [ i ]/h [ i +1]);
		
		B[0]=3*p[ 0 ] ;
	       	B[n-1]=3*p[ n-2]; //Gauss elimination :
		
		for ( int i =1;i<n ; i++){ D[ i] -= Q[ i -1]/D[ i -1]; B[ i] -= B[ i -1]/D[ i -1]; }
		
		b[ n-1]=B[ n-1]/D[ n-1]; //back−substitution :
		for ( int i=n-2;i >=0;i--) b[ i ]=(B[ i ]-Q[ i ]* b[ i +1])/D[ i ] ;
		for ( int i =0;i<n-1; i++){
			c[ i ]=(-2*b[ i ]-b[ i +1]+3*p[ i ] ) / h[ i ] ;
			d[ i ]=(b[ i ]+b[ i +1]-2*p[ i ] ) / h[ i ]/h[ i ] ;
		}	
	}

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
		int i = binsearch(x,z);

		double h=z-x[ i ]; // calculate the inerpolating spline :
		return y[ i ]+h*(b[ i ]+h*(c[ i ]+h*d[ i ] ) ) ;

	}
	
	public double derivative(double z){
		/* evaluate the derivative */
		//Derivative of si(x) is bi+2*ci*(x-xi)+3*di*(x-xi)^2

		int i = binsearch(x,z);
		double dsi = b[i]+2*c[i]*(z-x[i])+3*d[i]*Pow(z-x[i],2);
		return dsi;
	}//derivative

	public double integral(double z){
		/* evaluate the integral */
		int i = binsearch(x,z);

		double intOut = 0;
		for(int j=0;j<i;j++){  
		    intOut += (y[j]*(x[j+1]-x[j])+b[j]*Pow((x[j+1]-x[j]),2)/2+c[j]*Pow((x[j+1]-x[j]),3)/3)+d[j]*Pow((x[j+1]-x[j]),4)/4.0;
		}

		//Need to include point at z:
		intOut += (y[i]*(z-x[i])+b[i]*Pow((z-x[i]),2)/2+c[i]*Pow((z-x[i]),3)/3+d[i]*Pow((z-x[i]),4)/4);

		//integration gives us an additional constant which in this case is C = -1.0:
		double C = -1.0;
		intOut += C;

		return intOut;
		
	}//integral




public static int Main(){
		Random rnd = new Random(123);
		vector xs = new vector(11);	
		vector ys = new vector(11);

		//Func<double,double> test_function = x => Sin(2*x)*Exp(-x/3.0) + Cos(x);
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
			Console.Out.WriteLine($"{usex} {test_function(usex)}");
		}

		System.Console.Out.WriteLine("\n");
		///Cubic interp output
		var cspline = new cspline(xs,ys);
		
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {cspline.evaluate(usex)}");
		}

		System.Console.Out.WriteLine("\n");
		
		//Derivative output
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			System.Console.Out.WriteLine($"{usex} {cspline.derivative(usex)}");
		}
		
		System.Console.Out.WriteLine("\n");
	
		//Integral output
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			System.Console.Out.WriteLine($"{usex} {cspline.integral(usex)}");
		}
	
		System.Console.Out.WriteLine("\n");

		//We want to check that the built-in cubic splines in gnuplot produces a similar cubic spline, so we now try a more fun function than Sin:
		Func<double,double> fun_function = x => Sin(2*x)*Exp(-x/3.0)+Cos(x);

		for(int i=0;i<xs.size;i++){
			xs[i]=i;
			ys[i]=fun_function(xs[i]);
			Console.WriteLine($"{xs[i]} {ys[i]}");
		}
		
		System.Console.Out.WriteLine("\n");

		//Fun function output
		for(int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {fun_function(usex)}");
		}

		System.Console.Out.WriteLine("\n");
		///Cubic interp output for fun function
		var cspline_fun = new cspline(xs,ys);
		
		for (int i=0;i<numPoints;i++){
			double usex = start + i*step;
			Console.Out.WriteLine($"{usex} {cspline_fun.evaluate(usex)}");
		}

		System.Console.Out.WriteLine("\n");
		return 0;
	}
}

