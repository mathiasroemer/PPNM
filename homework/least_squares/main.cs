using System;
using static System.Math;

public class leastsq {
	public static (vector,matrix) lsfit(Func<double,double>[] fs, vector x, vector y, vector dy){
		int n = x.size, m=fs.Length;
		var A = new matrix(n,m);
		var b = new vector(n);
		for (int i =0;i<n ; i++){
			b[i] = y[i]/dy[i];
			for (int k=0;k<m;k++) A[i,k] = fs[k](x[i])/dy[i];
		}
		var QR = QRGS.decomp(A);
		var Q = QR.Item1;
		var R = QR.Item2;

		var c = QRGS.solve(Q,R,b);

		//S = (A^T * A)^-1, my inverse function requires Q and R as arguments, so A^T * A must first be in QR-decomposition such that the inverse can be found.

		var QRATA = QRGS.decomp(A.T * A);
		var QATA = QRATA.Item1;
		var RATA = QRATA.Item2;
		var ATAInv = QRGS.inverse(QATA,RATA);

		var S = ATAInv;
		return (c, S);
	}
	public static int Main() {
		System.Console.Out.WriteLine("Lets first be sure that the QR decomp works:");
		var rnd = new Random(123);
		
		var mat = matrix.random(8,3,rnd);
		mat.print("A (Random tall matrix):");
		
		var QR = QRGS.decomp(mat);
		var Q = QR.Item1;
		var R = QR.Item2;
		
		Q.print("Q:");
		R.print("R:");

		var QprodR = Q*R;
		QprodR.print("Q * R");

		System.Console.Out.WriteLine("We see R is upper triangular and the product of Q and R is A giving us that the QR decomposition routines work for tall matrices.\n");
		
		
		System.Console.Out.WriteLine("Lets now investigate the law of radioactive decay:\n");
		
		var fs = new Func<double,double>[] {z => 1, z => -z };
		vector x = new vector (new double[] {1,2,3,4,6,9,10,13,15});
		
		double[] ys = {117,100,88,72,53,29.5,25.2,15.2,11.1};
		double[] dys = {6,5,4,4,4,3,3,2,2};

		double[] lnys = new double[ys.Length];
		double[] lndys = new double[dys.Length];
		
		for(int i=0;i<dys.Length;i++){
			lnys[i]=Log(ys[i]);
			lndys[i]=Log(dys[i])/dys[i];
		}

		vector lny = new vector(lnys);
		vector lndy = new vector(lndys);

		var fit = lsfit(fs,x,lny,lndy);

		var c = fit.Item1;
		System.Console.Out.WriteLine("By doing an ordinary least squares fit on the data using ln(y) = ln(a)-l*t, we find the coefficients:");
		var a = Exp(c[0]);
		var lambda = c[1];
		System.Console.Out.WriteLine($"a = {a}, lambda = {lambda}\n");
		var T_12 = Log(2)/lambda;
		System.Console.Out.WriteLine($"Half-life time found: T_1/2 = Log(2)/lambda = {T_12} days");
		System.Console.Out.WriteLine("Half-life time modern value: T_1/2 = 3.631(2) days");
				
		var S = fit.Item2;
		S.print("The covariance matrix is found as:");

		var lambdaErr = Sqrt(S[1,1]);
		var T_12Err = lambdaErr * Log(2)/(Pow(lambda,2));
		System.Console.Out.WriteLine($"\nWith use of covariance matrix and error propagation, the error on T_1/2 is found as: {T_12Err}");

		System.Console.Out.WriteLine($"such that Half-life time found: T_1/2 = {T_12} +- {T_12Err} days");
		System.Console.Out.WriteLine("We see that the found half-life time lies within the uncertainty of the modern value.");
		return 0;	
	}
}
