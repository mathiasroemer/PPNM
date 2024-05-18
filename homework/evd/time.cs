using System;
using static System.Math;

public static class main {
	public static int Main(string[] args){
		int N = 1;
		foreach(var arg in args){
			var words = arg.Split(':');
			if(words[0]=="-N"){
				N = int.Parse(words[1]);
			}
		}
		//System.Console.WriteLine($"Make {N}x{N} matrix");
		
		var rnd = new Random(321);
		var A = matrix.randomSym(N,rnd);
		//System.Console.WriteLine($"Diagonalize {N}x{N} matrix");

		var eig = jacobi.cyclic(A);

		//System.Console.WriteLine($"Done {N}x{N} matrix");
		return 0;
	}
}
