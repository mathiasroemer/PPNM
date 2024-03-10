using System;
using static System.Math;

public static class jacobi {
	public static void timesJ(matrix A, int p, int q, double theta){
		double c=Cos(theta),s=Sin(theta);
		for(int i=0;i<A.size1;i++){
			double aip=A[i,p],aiq=A[i,q];
			A[i,p]=c*aip-s*aiq;
			A[i,q]=s*aip+c*aiq;
		}
	}
	public static void Jtimes(matrix A, int p, int q, double theta){
		double c=Cos(theta),s=Sin(theta);
		for(int j=0;j<A.size1;j++){
			double apj=A[p,j],aqj=A[q,j];
			A[p,j]= c*apj+s*aqj;
			A[q,j]=-s*apj+c*aqj;
		}
	}
	public static (vector,matrix) cyclic(matrix M){
		matrix A=M.copy();
		matrix V=matrix.id(M.size1);
		vector w=new vector(M.size1);
		
		int n = A.size1;

		/* run Jacobi rotations on A and update V */
		bool changed;
		do{
			changed=false;
			for(int p=0;p<n-1;p++)
			for(int q=p+1;q<n;q++){
				double apq=A[p,q], app=A[p,p], aqq=A[q,q];
				double theta=0.5*Atan2(2*apq,aqq-app);
				double c=Cos(theta),s=Sin(theta);
				double new_app=c*c*app-2*s*c*apq+s*s*aqq;
				double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
				if(new_app!=app || new_aqq!=aqq) // do rotation
					{
					changed=true;
					timesJ(A,p,q, theta); // A←A*J 
					Jtimes(A,p,q,-theta); // A←JT*A 
					timesJ(V,p,q, theta); // V←V*J
					}
			}
		}while(changed);
		
		/* copy diagonal elements into w */
		for(int i=0;i<n;i++){
			w[i]=A[i,i];
		}
		return (w,V);
	}
	
	public static int Main(){
		var rnd = new Random(321);
		var A = matrix.randomSym(4,rnd);
		A.print("A (Random symmetric matrix):");
		System.Console.WriteLine("");

		var eig = cyclic(A);
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
		
		return 0;
	}
}
