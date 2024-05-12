using System;
using static System.Math;

public static class main {
	public static int Main(string[] args){
		double rmax = 0;
		double dr = 0;
		bool egfnc = false;

		foreach(var arg in args){
			var words = arg.Split(':');
			if(words[0]=="-rmax"){
				rmax = double.Parse(words[1]);
			}
			else if(words[0]=="-dr"){
				dr = double.Parse(words[1]);
			}
			else if(words[0]=="-egfnc"){
				if(words[1]=="y") egfnc = true;
			}
		}
		
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
		
		var H_eig = jacobi.cyclic(H);
		var H_w = H_eig.Item1;
		var H_V = H_eig.Item2;
		
		double eps_0 = H_w[0];
		
		if(egfnc == false) {
			System.Console.WriteLine($"{rmax} {dr} {eps_0}");
		}
		else {
			for(int i=0;i<2;i++){
				for(int j=0;j<npoints;j++) {
					System.Console.WriteLine($"{r[j]} {H_V[j,i]/Sqrt(dr)}");
				}
			System.Console.WriteLine("\n");
		}
		}
		return 0;
	}
}
