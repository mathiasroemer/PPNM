using System;
using static System.Math;

public class minimization{
	public static vector gradient(Func<vector,double> φ,vector x){
		vector gφ = new vector(x.size);
		double φx = φ(x); /* no need to recalculate at each step */
		for(int i=0;i<x.size;i++){
			double dx=Max(Abs(x[i]),1)*Pow(2,-26);
			x[i]+=dx;
			gφ[i]=(φ(x)-φx)/dx;
			x[i]-=dx;
		}
		return gφ;
	}

	public static matrix hessian(Func<vector,double> φ,vector x){
		matrix H=new matrix(x.size);
		vector gφx=gradient(φ,x);
		for(int j=0;j<x.size;j++){
			double dx=Max(Abs(x[j]),1)*Pow(2,-26); /* for numerical gradient */
			x[j]+=dx;
			vector dgφ=gradient(φ,x)-gφx;
			for(int i=0;i<x.size;i++) H[j,i]=dgφ[i]/dx;
			x[j]-=dx;
		}
		return H;
		//return (H+H.T)/2; // you think?
	}

	public static (vector,int) newton(
		Func<vector,double> φ, /* objective function */
		vector x,              /* starting point */
		double acc=1e-5        /* accuracy goal, on exit |∇φ| should be < acc */
	){
		int NSteps = 0;
		int maxsteps = 1000;
		double λmin = 1/64.0;
		do{ /* Newton's iterations */
			NSteps++;
			var gφ = gradient(φ,x);
			if(gφ.norm() < acc) break; /* job done */
			var H = hessian(φ,x);
			//var QRH = givensQR(H);   /* QR decomposition */
			//var dx = QRH.solve(-∇φ); /* Newton's step */
			var QRH = QRGS.decomp(H);
			var (Q,R) = QRH;
			vector dx = QRGS.solve(Q,R,-gφ);
		
			double λ=1.0,φx=φ(x);
			//double λmin = 1.0/64.0;
			
			do{ /* linesearch */
				if( φ(x+λ*dx) < φx ) break; /* good step: accept */
				if( λ < λmin ) break; /* accept anyway */
				λ/=2.0;
			}while(true);
			x+=λ*dx;
		}while(NSteps < maxsteps);
		return (x,NSteps);
	}//newton
}
