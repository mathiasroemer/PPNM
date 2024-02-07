using System;
using static vec;
using static System.Console;
using static System.Math;

class main{
        static void Main(){
		var a = new vec(0,2,3);
		var b = new vec(3,1,0);
		var c = new vec(2,0,1);
		vec[] vecs = {a,b,c};
		int i = 1;
		WriteLine("Vectors:");
		foreach (vec cVec in vecs){
			cVec.print("Vec "+i+": ");
			i++;
		}
		WriteLine("\nVec 1 . Vec 2, Vec 2 . Vec 3 and Vec 3 . Vec 1:");
		WriteLine($"{vec.dot(a,b)}, {vec.dot(b,c)} and {vec.dot(c,a)}\n");
		WriteLine("Vec 1 + 2, Vec 2 + 3 and Vec 3 + 1:");
		WriteLine($"{a+b}, {b+c} and {c+a}\n");
		double cnst = 3.0;
		WriteLine($"constant = {cnst}, c * Vec 1, c * Vec 2 and c * Vec 3:");
		WriteLine($"{cnst*a}, {cnst*b} and {cnst*c}\n");
		WriteLine("Vec 1 approx 2, vec 2 approx 3 and vec 3 approx 1:");
		WriteLine($"{a.approx(b)}, {b.approx(c)} and {c.approx(a)}");
	}
}
