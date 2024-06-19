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
}

