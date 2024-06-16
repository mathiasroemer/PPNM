using System;
using static System.Math;

public class roots{
	public static matrix jacobian
	(Func<vector,vector> f,vector x,vector fx=null,vector dx=null){
		if(dx == null) dx = x.map(xi => Abs(xi)*Pow(2,-26));
		if(fx == null)fx=f(x);
		matrix J=new matrix(x.size);
		for(int j=0;j < x.size;j++){
			x[j]+=dx[j];
			vector df=f(x)-fx;
			for(int i=0;i < x.size;i++) J[i,j]=df[i]/dx[j];
			x[j]-=dx[j];
			}
		return J;
	}//jacobian
	
	public static vector newton(
	Func<vector,vector>f /* the function to find the root of */
	,vector x            /* the start point */
	,double acc=1e-2     /* accuracy goal: on exit ‖f(x)‖ should be <acc */
	,vector δx=null      /* optional δx-vector for calculation of jacobian */
	){
		vector fx=f(x),z,fz;
		do{ /* Newton's iterations */
			if(fx.norm() < acc) break; /* job done */
			matrix J=jacobian(f,x,fx,δx);
			var QRofJ = QRGS.decomp(J);
			var (Q,R) = QRofJ;
			vector Dx = QRGS.solve(Q,R,-fx);
			double λ=1;
			double λmin=1.0/64; 
			do{ /* linesearch */
				z=x+λ*Dx;
				fz=f(z);
				if( fz.norm() < (1-λ/2)*fx.norm() ) break;
				if( λ < λmin ) break;
				λ/=2;
				}while(true);
			x=z; fx=fz;
			}while(true);
		return x;
	}//newton

	public static vector newtonQuad(
		Func<vector,vector>f /* the function to find the root of */
		,vector x            /* the start point */
		,double acc=1e-2     /* accuracy goal: on exit ‖f(x)‖ should be <acc */
		,vector δx=null      /* optional δx-vector for calculation of jacobian */
		){
			vector fx=f(x),z,fz;
			do{ /* Newton's iterations */
				if(fx.norm() < acc) break; /* job done */
				matrix J=jacobian(f,x,fx,δx);
				var QRofJ = QRGS.decomp(J);
				var (Q,R) = QRofJ;
				vector Dx = QRGS.solve(Q,R,-fx);
				double λ=1.0;
				double λmin=1.0/64; 
				do{ /* quadratic linesearch */
					z=x+λ*Dx;
					fz=f(z);
					
					if(fz.norm() < (1-λ/2)*fx.norm())  break;

					if( λ < λmin ) break;
						
					double phizero = (1.0/2.0) * Math.Pow((fx.norm()),2);
					double phidzero = -Math.Pow(fx.norm(),2);

					double c = ((1.0/2.0) * Math.Pow(fz.norm(),2) - phizero - phidzero*λ)/Math.Pow(λ,2); 

					λ = -phidzero/(2.0*c);
					}while(true);
				x=z; fx=fz;
				}while(true);
			return x;
	}//newtonQuad
}
