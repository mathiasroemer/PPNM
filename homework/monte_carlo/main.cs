using System;
using static System.Math;
using System.Diagnostics;

public class main{
	public static (double,double) plainmc(Func<vector,double> f,vector a,vector b,int N){
		int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
		double sum=0,sum2=0;
		var x=new vector(dim);
		var rnd=new Random();
		for(int i=0;i<N;i++){
			for(int k=0;k<dim;k++)x[k]=a[k]+rnd.NextDouble()*(b[k]-a[k]);
			double fx=f(x); sum+=fx; sum2+=fx*fx;
			}
		double mean=sum/N, sigma=Sqrt(sum2/N-mean*mean);
		var result=(mean*V,sigma*V/Sqrt(N));
		return result;
	}

	public static double corput(int n, int b){
		double q=0, bk=(double)1/b ;
		while(n>0){ q += (n % b)*bk ; n /= b ; bk /= b ; }
		return q ; 
	}
	
	public static void halton(int n, vector x){
		int[] bases=new int[] {2 ,3 ,5 ,7 ,11 ,13 ,17 ,19 ,23 ,29 ,31 ,37 ,41 ,43 ,47 ,53 ,59 ,61};
		int maxd = bases.Length;
		int d = x.size;
		Debug.Assert (d <= maxd); 
		for ( int i =0;i<d ; i++) { x [ i ]=corput (n , bases [ i ] ) ; }
	
	}

	public static void lattice(int n, vector x){
		int d = x.size;
		vector alpha = new vector(d);
		for (int i=0;i<d;i++){
			alpha[i]=Sqrt((Math.PI+i)) % 1;
		}

		for (int i=0;i<d;i++) {
			x[i] = (n*alpha[i]) % 1;
		}
	}

	public static (double,double) quasimc(Func<vector,double> f,vector a,vector b,int N,int type){
		int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
		double sum=0,sum2=0;
		var x=new vector(dim);
		var haltonVec = new vector(dim);
		var latticeVec = new vector(dim);

		for(int i=0;i<N;i++){
			halton(i,haltonVec);
			lattice(i,latticeVec);
			for(int k=0;k<dim;k++)x[k]=a[k]+haltonVec[k]*(b[k]-a[k]);
			sum+=f(x);
			
			for(int k=0;k<dim;k++)x[k]=a[k]+latticeVec[k]*(b[k]-a[k]);
			sum2+=f(x);
		}
		double mean=0;
		if (type==1){
			mean=sum/N;
		}
		else if (type==2){
			mean=sum2/N;
		}
		double sigma=Abs(sum-sum2)/N;
		var result=(mean*V,sigma*V);
		return result;
	}

	public static int Main(){
		Console.WriteLine("Part A:\n");
		
		double R = 1.0;
		vector start1 = new vector(0,0);
		vector end1 = new vector(R,2*Math.PI);
		
		Console.WriteLine("To calculate the area of a circle in polar coordinates, we let r go from 0 to R, and theta from 0 to 2*PI, and then we just integrate over r.");
		Func<vector,double> func = delegate(vector v){return v[0];};
	
		int testN = 10000;
		var result1 = plainmc(func,start1,end1,testN);
		Console.WriteLine($"Using the plain Monte Carlo Integration we get an area of a unit circle using {testN} points: \n{result1.Item1:F5} +- {result1.Item2:F5}");

		for(int N=10;N<2.5e4;N+=500) {
			var res = plainmc(func,start1,end1,N);
			Console.Error.WriteLine($"{N} {res.Item1} {res.Item2} {Abs(res.Item1 - Math.PI)}");
		}
		
		Console.WriteLine("\n");

		Console.WriteLine("We now try to calculate the difficult singular integral using 1000000 points:");

		vector start2 = new vector(0,0,0);
		vector end2 = new vector(Math.PI,Math.PI,Math.PI);
		Func<vector,double> difficultFunc = delegate(vector v){ return (Math.Pow(1-(Math.Cos(v[0])*Math.Cos(v[1])*Math.Cos(v[2])),-1))/Math.Pow(Math.PI,3);};

		var result2 = plainmc(difficultFunc,start2,end2,1000000);
		Console.WriteLine($"{result2.Item1:F5} +- {result2.Item2:F5}");
		
		Console.WriteLine("\nPart B:\n");
		Console.WriteLine("We now introduce the Van der Corput, Halton and lattice rule sequences.");
		Console.WriteLine("We first use the Halton sequence to generate quasi-random sampling for the integration, and then use the lattice rule sequence to get an error estimate by using the difference of the two sequences.");
		
		var result3 = quasimc(func,start1,end1,testN,1);
		Console.WriteLine($"Using the Quasi-random Monte Carlo Integration, with the Halton sequence, we get an area of a unit circle using {testN} points: \n{result3.Item1:F5} +- {result3.Item2:F5}");

		var result4 = quasimc(func,start1,end1,testN,2);
		Console.WriteLine($"Using the Quasi-random Monte Carlo Integration, with the lattice rule sequence, we get an area of a unit circle using {testN} points: \n{result4.Item1:F5} +- {result4.Item2:F5}");

		Console.Error.WriteLine("\n");
		for(int N=10;N<2.5e4;N+=500) {
			var res = quasimc(func,start1,end1,N,1);
			Console.Error.WriteLine($"{N} {res.Item1} {res.Item2}");
		}

		Console.Error.WriteLine("\n");
		for(int N=10;N<2.5e4;N+=500) {
			var res1 = plainmc(func,start1,end1,N);
			var res2 = quasimc(func,start1,end1,N,1);
			Console.Error.WriteLine($"{N} {Abs(res1.Item1 - Math.PI)} {Abs(res2.Item1 - Math.PI)}");
		}

		return 0;
	}
}
