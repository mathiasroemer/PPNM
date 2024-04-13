using System;
using static System.Math;

public class main{
	public static (vector,vector) rkstep12(
		Func<double,vector,vector> f,/* the f from dy/dx=f(x,y) */
		double x,                    /* the current value of the variable */
		vector y,                    /* the current value y(x) of the sought function */
		double h                     /* the step to be taken */
		)
	{
		vector k0 = f(x,y);              /* embedded lower order formula (Euler) */
		vector k1 = f(x+h/2,y+k0*(h/2)); /* higher order formula (midpoint) */
		vector yh = y+k1*h;              /* y(x+h) estimate */
		vector δy = (k1-k0)*h;           /* error estimate */
		return (yh,δy);
	}

	public static (genlist<double>,genlist<vector>) driver(
	Func<double,vector,vector> F,/* the f from dy/dx=f(x,y) */
	(double,double) interval,    /* (start-point,end-point) */
	vector ystart,               /* y(start-point) */
	double h=0.125,              /* initial step-size */
	double acc=0.01,             /* absolute accuracy goal */
	double eps=0.01              /* relative accuracy goal */
	){
	var (a,b)=interval; double x=a; vector y=ystart.copy();
	var xlist=new genlist<double>(); xlist.add(x);
	var ylist=new genlist<vector>(); ylist.add(y);
	do{
		if(x>=b) return (xlist,ylist); /* job done */
		if(x+h>b) h=b-x;               /* last step should end at b */
		var (yh,δy) = rkstep12(F,x,y,h);
		double tol = (acc+eps*yh.norm()) * Sqrt(h/(b-a));
		double err = δy.norm();
		if(err<=tol){ // accept step
			x+=h; y=yh;
			xlist.add(x);
			ylist.add(y);
			}
		h *= Min( Pow(tol/err,0.25)*0.95 , 2); // readjust stepsize
		}while(true);
	}//driver

	public static int binsearch(double[] x, double z)
		{/* locates the interval for z by bisection */ 
		if( z<x[0] || z>x[x.Length-1] ) throw new Exception("binsearch: bad z");
		int i=0, j=x.Length-1;
		while(j-i>1){
			int mid=(i+j)/2;
			if(z>x[mid]) i=mid; else j=mid;
			}
		return i;
	}//binsearch

	public static Func<double,vector> make_linear_interpolant(genlist<double> x,genlist<vector> y)
	{
		Func<double,vector> interpolant = delegate(double z){
			int i=binsearch(x,z);
			double Δx=x[i+1]-x[i];
			vector Δy=y[i+1]-y[i];
			return y[i]+Δy/Δx*(z-x[i]);
		};
		return interpolant;
	}//make_linear_interpolant

	public static Func<double,vector> make_ode_ivp_interpolant
	(Func<double,vector,vector> f,(double,double)interval,vector y,double acc=0.01,double eps=0.01,double hstart=0.01 )
	{
		var (xlist,ylist) = driver(f,interval,y,acc,eps,hstart);
		return make_linear_interpolant(xlist,ylist);
	}//make_ode_ivp_interpolant

	public static int Main(){
		System.Console.WriteLine("Part A:\n");

		System.Console.WriteLine("We first test the routines on u''=-u\n");
		Func<double, vector, vector> F = delegate(double x, vector y){
			return new vector(y[1],-y[0]);
		};
		
		vector ystart = new vector(0,1);

		var (xlist,ylist) = driver(F,(0,10.0),ystart);

		for(int i=0;i<xlist.size;i++) {
			System.Console.Error.WriteLine($"{xlist[i]} {ylist[i][0]} {ylist[i][1]}");
		}

		System.Console.Error.WriteLine("\n");

		System.Console.WriteLine("We then try to solve the Lotka-Volterra system\n");
		double a = 1.5; double b = 1.0; double c = 3.0; double d = 1.0;
		Func<double, vector, vector> F_LV = delegate(double t, vector z) {
			var x = z[0]; var y = z[1];
			return new vector(a*x-b*x*y,-c*y+d*x*y);
		};

		vector yinit = new vector(10,5);
		var (xs,ys) = driver(F_LV,(0,10),yinit);
		for(int i=0;i<xs.size;i++){
			System.Console.Error.WriteLine($"{xs[i]} {ys[i][0]} {ys[i][1]}");
		}
		System.Console.Error.WriteLine("\n");

		System.Console.WriteLine("Part B:\n");

		System.Console.WriteLine("We now want to use the routines on u'' + u= l + epsilon * u^2\n");
		
		double epsilon = 0;
		int NPoints = 1000;
		double EndPoint = 15*System.Math.PI;

		Func<double, vector, vector> F_eqmot = delegate(double x, vector y){
			return new vector(y[1],1-y[0]+epsilon*y[0]*y[0]);
		};
			
		vector ystart1 = new vector(1,0);

		var interp1 = make_ode_ivp_interpolant(F_eqmot,(0,EndPoint),ystart1);
		for(int i = 0; i<NPoints; i++){
			double x = (i*EndPoint)/NPoints;
			System.Console.Error.WriteLine($"{x} {interp1(x)[0]} {interp1(x)[1]}");	
		}
		
		System.Console.Error.WriteLine("\n");

		vector ystart2 = new vector(1,-0.5);
		epsilon = 0;

		var interp2 = make_ode_ivp_interpolant(F_eqmot,(0,EndPoint),ystart2);
		for(int i = 0; i<NPoints; i++){
			double x = (i*EndPoint)/NPoints;
			System.Console.Error.WriteLine($"{x} {interp2(x)[0]} {interp2(x)[1]}");	
		}
		
		System.Console.Error.WriteLine("\n");

		vector ystart3 = new vector(1,-0.5);
		epsilon = 0.01;

		var interp3 = make_ode_ivp_interpolant(F_eqmot,(0,EndPoint),ystart3);
		for(int i = 0; i<NPoints; i++){
			double x = (i*EndPoint)/NPoints;
			System.Console.Error.WriteLine($"{x} {interp3(x)[0]} {interp3(x)[1]}");	
		}
		
		return 0;
	}
}
