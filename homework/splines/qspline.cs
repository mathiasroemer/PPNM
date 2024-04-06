using System;
using static System.Math;

public class qspline{
	vector x,y,b,c;
	public qspline(vector xs, vector ys){
		/* copy xs,ys into x,y and build b,c */
		this.x = xs.copy();
		this.y = ys.copy();
		System.Console.Out.WriteLine("qspline init");
	}
	public double evaluate(double z){
		/* evaluate the spline using x,y,b,c */
		return 1.0;
	}
	public static int Main(){
		System.Console.Out.WriteLine("asd");
		return 0;
	}
}


