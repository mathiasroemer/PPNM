using System;
using static System.Math;

public class bettersubspline {
	vector x,y,dy,a,b,c,d,e;
	public bettersubspline(vector xs,vector ys,vector dys){
		x = xs.copy();
		y = ys.copy();
		dy = dys.copy();
	
		a = new vector(x.size);
		b = new vector(x.size);
		c = new vector(x.size-1);
		d = new vector(x.size-1);
		e = new vector(x.size-1);

		Func<int,double> diffxi = i => x[i+1]-x[i];
		Func<int,double> diffyi = i => y[i+1]-y[i];
		Func<int,double> dyi = i => dy[i];
		Func<int,double> pi = i => diffyi(i)/diffxi(i);

		for(int i=0;i<a.size;i++) a[i] = y[i];
		for(int i=0;i<b.size;i++) b[i] = dyi(i);
		for(int i=0;i<c.size;i++) c[i] = (3*pi(i)-2*dyi(i)-dyi(i+1))/(diffxi(i));
		for(int i=0;i<d.size;i++) d[i] = (dyi(i)+dyi(i+1)-2*pi(i))/(Math.Pow(diffxi(i),2));
		
		//Forward-recursion:
		e[0] = 0;

		for(int i=0;i<e.size-1;i++){
			e[i+1] = (1.0/(Math.Pow(diffxi(i+1),2))) * (c[i]-c[i+1]+3*d[i]*diffxi(i)+e[i]*Math.Pow(diffxi(i),2));
		}

		//Backward-recursion:
		e[e.size-1] /= 2.0;
		for(int i=e.size-2;i>1;i--){
			e[i] = (1.0/(Math.Pow(diffxi(i),2))) * (c[i+1]-c[i]-3*d[i]*diffxi(i)+e[i+1]*Math.Pow(diffxi(i+1),2));
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
		double h=z-x[i]; // calculate the interpolating spline:
		return y[i]+h*(b[i]+h*(c[i]+h*d[i]))+e[i]*Math.Pow(z-x[i],2)*Math.Pow(z-x[i+1],2);

	}
	
	public double derivative(double z){
		/* evaluate the derivative */
		//Derivative of si(x) is bi+2*ci*(x-xi)+3*di*(x-xi)^2

		int i = binsearch(x,z);
		double dsi = b[i]+2*c[i]*(z-x[i])+3*d[i]*Pow(z-x[i],2)+2*e[i]*(z-x[i])*(z-x[i+1])*(2*z-x[i]-x[i+1]);
		return dsi;
	}//derivative
			
	public double second_derivative(double z){
		/* evaluate the second derivative */

		int i = binsearch(x,z);
		double ddsi = 2*c[i]+6*d[i]*(z-x[i])+2*e[i]*(6*Math.Pow(z,2)-6*z*(x[i]+x[i+1])+Math.Pow(x[i],2)+Math.Pow(x[i+1],2)+4*x[i]+x[i+1]);
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

