using System;
using static System.Math;

public static class main {
	public static int Main(string[] args){
		System.Console.WriteLine("Part A:\n");
		var rnd = new Random(321);
		var A = matrix.randomSym(4,rnd);
		A.print("A (Random symmetric matrix):");
		System.Console.WriteLine("");

		var eig = jacobi.cyclic(A);
		var w = eig.Item1;
		var V = eig.Item2;
		
		System.Console.WriteLine("Implement Jacobi eigenvalue algorithm on A:");
		w.print("w:");
		System.Console.WriteLine("");
		V.print("V:");
		System.Console.WriteLine("");
		
		var D = new matrix(w);
		D.print("We get D from w:");	
		System.Console.WriteLine("");
		
		var VTAV = V.T * A * V;
		VTAV.print("V^T * A * V:");
		System.Console.Write("Using matrix approx function with D we get: ");
		System.Console.WriteLine(VTAV.approx(D));
		System.Console.WriteLine("We see that it is approximately equal to D.");
		System.Console.WriteLine("");

		var VDVT = V * D * V.T;
		VDVT.print("V * D * V.T:");
		System.Console.WriteLine("We see that it is equal to A.\n");
		
		var VTV = V.T * V;
		VTV.print("V^T * V:");
		System.Console.WriteLine("We see that it is equal to the identity matrix.\n");
	
		var VVT = V * V.T;
		VVT.print("V * V.T:");
		System.Console.WriteLine("We see that it is also equal to the identity matrix.\n");
		
		System.Console.WriteLine("Part B:\n");
		double rmax = 0;
		double dr = 0;
		foreach(var arg in args){
			var words = arg.Split(':');
			if(words[0]=="-rmax"){
				rmax = double.Parse(words[1]);
			}
			else if(words[0]=="-dr"){
				dr = double.Parse(words[1]);
			}
		}
		System.Console.Write("Chosen rmax: ");
		System.Console.WriteLine(rmax);
		System.Console.Write("Chosen dr: ");
		System.Console.WriteLine(dr);
		
		System.Console.WriteLine("\nWe first build the Hamiltonian matrix");

		int npoints = (int)(rmax/dr)-1;
		vector r = new vector(npoints);
		for(int i=0;i<npoints;i++)r[i]=dr*(i+1);
		matrix H = new matrix(npoints,npoints);
		for(int i=0;i<npoints-1;i++){
		   H[i,i]  =-2*(-0.5/dr/dr);
		   H[i,i+1]= 1*(-0.5/dr/dr);
		   H[i+1,i]= 1*(-0.5/dr/dr);
		  }
		H[npoints-1,npoints-1]=-2*(-0.5/dr/dr);
		for(int i=0;i<npoints;i++)H[i,i]+=-1/r[i];
		
		System.Console.WriteLine("Then we use our Jacobi routine to obtain eigenvalues and eigenvectors\n");
		var H_eig = jacobi.cyclic(H);
		var H_w = H_eig.Item1;
		var H_V = H_eig.Item2;
		
		H_w.print("Eigenvalues of Hamiltonian matrix:");

		double eps_0 = H_w[0];
		System.Console.WriteLine($"We know that the eigenvalue diagonal matrix will be arranged in accending order meaning that index 0,0 will be the lowest eigenvalue for the Hamiltonian matrix with max radius {rmax:f2} and step {dr:f2} Bohr radii: {eps_0:f3} Hartree\n");
		
		System.Console.WriteLine("We can now try to vary rmax and dr investigating the convergence of the lowest eigenvalue\n");
		System.Console.WriteLine($"By varying rmax and dr we can also plot the 2 lowest eigen-functions corresponding to that rmax and dr and compare them to the analytic reduced radial wave function\n");
		
		System.Console.WriteLine("Part C:\n");
		System.Console.WriteLine("For part c we choose to look at how the number of operations for matrix diagonalization scales as O(n^3). The data is generated in parallel utilizing shell's background operator &. Afterwards the data is fitted to a function: f(x) = a*x^3 +b");
		return 0;
	}
}
