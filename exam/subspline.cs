using System;
using static System.Math;

public class subspline {
	vector x,y,dy,a,b,c,d;
	public subspline(vector xs,vector ys,vector dys){
		x = xs.copy();
		y = ys.copy();
		dy = dys.copy();
	
		a = new vector(x.size);
		b = new vector(x.size);
		c = new vector(x.size-1);
		d = new vector(x.size-1);

		Func<int,double> diffxi = i => x[i+1]-x[i];
		Func<int,double> diffyi = i => y[i+1]-y[i];
		Func<int,double> dyi = i => dy[i];
		Func<int,double> pi = i => diffyi(i)/diffxi(i);

		for(int i=0;i<a.size;i++) a[i] = y[i];
		for(int i=0;i<b.size;i++) b[i] = dyi(i);
		for(int i=0;i<c.size;i++) c[i] = (3*pi(i)-2*dyi(i)-dyi(i+1))/(diffxi(i));
		for(int i=0;i<d.size;i++) d[i] = (dyi(i)+dyi(i+1)-2*pi(i))/(Math.Pow(diffxi(i),2));
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

		double h=z-x[i]; // calculate the interpolating spline:
		return y[i]+h*(b[i]+h*(c[i]+h*d[i]));

	}
	
	public double derivative(double z){
		/* evaluate the derivative */
		//Derivative of si(x) is bi+2*ci*(x-xi)+3*di*(x-xi)^2

		int i = binsearch(x,z);
		double dsi = b[i]+2*c[i]*(z-x[i])+3*d[i]*Pow(z-x[i],2);
		return dsi;
	}//derivative
			
	public double second_derivative(double z){
		/* evaluate the second derivative */

		int i = binsearch(x,z);
		double ddsi = 2*c[i]+6*d[i]*(z-x[i]);
		return ddsi;
	}//second derivative


	public double integral(double z){
		/* evaluate the integral */
		int i = binsearch(x,z);

		double intOut = 0;
		for(int j=0;j<i;j++){  
		    intOut += (y[j]*(x[j+1]-x[j])+b[j]*Pow((x[j+1]-x[j]),2)/2+c[j]*Pow((x[j+1]-x[j]),3)/3)+d[j]*Pow((x[j+1]-x[j]),4)/4.0;
		}

		//Need to include point at z:
		intOut += (y[i]*(z-x[i])+b[i]*Pow((z-x[i]),2)/2+c[i]*Pow((z-x[i]),3)/3+d[i]*Pow((z-x[i]),4)/4);

		//integration gives us an additional constant which in this case is C = 0:
		double C = 0;
		intOut += C;


		return intOut;
		
	}//integral
}

