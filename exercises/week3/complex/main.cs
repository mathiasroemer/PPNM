using System;
using static System.Console;

class main {
	static void Main(){
		WriteLine("sqrt(-1):");
		complex comp = new complex(-1,0);
		var sqcomp = cmath.sqrt(comp);
		cmath.print(sqcomp,"sqrt(-1) = ");
		cmath.print(cmath.I,"I = ");
		WriteLine($"sqrt(-1) approx I : {sqcomp.approx(cmath.I)}");
		WriteLine($"sqrt(-1) approx -I : {sqcomp.approx(-cmath.I)}");

		WriteLine("");
		WriteLine("ln(i):");
		var lni = cmath.log(cmath.I);
		cmath.print(lni,"ln(i) = ");
		var ipi2 = cmath.I*System.Math.PI/2;
		cmath.print(ipi2,"I*PI/2 = ");
		WriteLine($"ln(i) approx I*PI/2 : {lni.approx((cmath.I*System.Math.PI/2))}");
		
		WriteLine("");
		WriteLine("sqrt(i):");
		var sqrti = cmath.sqrt(cmath.I);
		cmath.print(sqrti,"sqrt(i) = ");
		var oosq2 = new complex(1/cmath.sqrt(2),1/cmath.sqrt(2));
		cmath.print(oosq2,"1/sqrt(2) + i/sqrt(2) = ");
		WriteLine($"sqrt(i) approx 1/sqrt(2) + i/sqrt(2): {sqrti.approx(oosq2)}");

		WriteLine("");
		WriteLine("i^i:");
		var ii = cmath.pow(cmath.I,cmath.I);
		cmath.print(ii,"i^i = ");
		var real = 0.208;
		WriteLine($"sqrt(i) approx 0.208: {ii.approx(real)}");
	}
}
