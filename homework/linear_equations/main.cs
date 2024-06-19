using System;

public static class QRGS{
   public static (matrix,matrix) decomp(matrix A){
	matrix Q=A.copy();
	int m = A.size2;
	matrix R=new matrix(m,m);
      	/* orthogonalize Q and fill-in R */
	for(int i=0;i<m;i++){
		R[i,i] = Q[i].norm();
		Q[i] /= R[i,i];
		for(int j=i+1;j<m;j++){
			R[i,j]=Q[i].dot(Q[j]);
			Q[j] -= Q[i]*R[i,j];
		}
	}
      return (Q,R);
      }
   public static vector solve(matrix Q, matrix R, vector b){
	//using R*x = Q^T * b and backsubstituion:
	var c = Q.T * b;
	var U = R.copy();
	for(int i=c.size-1;i >= 0; i--){
		double sum=0;
		for(int k=i+1;k<c.size;k++) sum+=U[i,k]*c[k];
		c[i] = (c[i]-sum)/U[i,i];
	}
	return c;	
   }
   public static double det(matrix R){
   	double det = 1;
	int m = R.size2;

	for (int i = 0; i<m;i++) det *= R[i,i];
	return det;
   }
   public static matrix inverse(matrix Q,matrix R){
	   int n = Q.size1;
	   
	   matrix Ainv = new matrix(n,n);
	   vector unit = new vector(n);

	   for(int i=0;i<n;i++){
		   unit[i] = 1;
		   Ainv[i] = solve(Q,R,unit);
		   unit[i] = 0;
	   }
	   
	   return Ainv;
   }
   public static int Main(){
	var rnd = new Random(321);

	var A = matrix.random(6,3,rnd);
	A.print("A (Random matrix):");
	
	System.Console.WriteLine("\nLets check decomp:\n");

	var QRfact = decomp(A);
	var Q = QRfact.Item1;
	var R = QRfact.Item2;
	
	Q.print("Q:");
	System.Console.WriteLine("");

	R.print("R:");
	System.Console.WriteLine("R is upper triangular\n");

	var QTQ = Q.T * Q;
	QTQ.print("Q^T * Q:");
	System.Console.WriteLine("Q^T * Q seems to be approximately equal to 1 (the identity matrix)\n");

	var QR = Q * R;
	QR.print("Q * R:");
	System.Console.WriteLine("We see that Q * R is equal to A\n");

	System.Console.WriteLine("\nLets check solve:\n");
	
	var squareMat = matrix.random(5,5,rnd);
	var b = vector.random(5,rnd);

	squareMat.print("A (Random square matrix):");
	System.Console.WriteLine("");
	
	b.print("b (Random vector):");
	System.Console.WriteLine("");

	var newQRfact = decomp(squareMat);
	var newQ = newQRfact.Item1;
	var newR = newQRfact.Item2;
	
	newQ.print("Q:");
	System.Console.WriteLine("");
	
	newR.print("R:");
	System.Console.WriteLine("R is again upper triangular");
	System.Console.WriteLine("");

	var x = solve(newQ,newR,b);
	x.print("Solving QRx=b we get, x:");
	System.Console.WriteLine("");

	var Ax = squareMat * x;
	Ax.print("A*x:");
	System.Console.WriteLine("Which is equal to b\n");

	System.Console.WriteLine("Lets find the inverse of the random square matrix A:");
	var B = inverse(newQ,newR);
	B.print("B (A inverse):");
	System.Console.WriteLine("");

	var AB = squareMat*B;
	AB.print("A*B:");
	System.Console.WriteLine("We see that A*B is approximately the identity matrix");
	return 0;
   }
}
