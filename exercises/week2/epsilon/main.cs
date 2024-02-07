using System;
using static System.Console;
using static System.Math;

class epsilon {
	public static bool approx(double a, double b, double acc=1e-9, double eps=1e-9){
		if(Abs(b-a) <= acc) return true;
		if(Abs(b-a) <= Max(Abs(a),Abs(b))*eps) return true;
		return false;
	}
	static void Main(){
		int i=1; while(i+1>i) {i++;}
		Write("my max int = {0}\n",i);
		Write($"int.MaxValue = {int.MaxValue}\n");
		
		i=1; while(i-1<i) {i++;}
		Write("my min int = {0}\n",i);
		Write($"int.MinValue = {int.MinValue}\n");
		double x=1; while(1+x!=1){x/=2;} x*=2;
		float y=1F; while((float)(1F+y) != 1F){y/=2F;} y*=2F;
		WriteLine("\nMachine epsilon for double and float:");
		WriteLine($"double x = {x}, Pow(2,-52) = {Pow(2,-52)}");
		WriteLine($"float y = {y}, Pow(2,-23) = {Pow(2,-23)}");
		WriteLine("\nTiny:");
		double epsilon=Pow(2,-52);
		double tiny=epsilon/2;
		double a=1+tiny+tiny;
		double b=tiny+tiny+1;
		Write($"a==b ? {a==b} - They are not equal since the order of operations matter in this case.\n");
		Write($"a>1  ? {a>1} - This is false since the tiny won't be enough to make a precision rounding to above the 1.\n");
		Write($"b>1  ? {b>1} - This is true since the two tiny's will round such that they are greater than 0, we than add 1 such that b will be greater than 1.\n");
		double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
		double d2 = 8*0.1; 
		WriteLine("\n0.1 and == operator:");
		WriteLine($"d1={d1:e15}");
		WriteLine($"d2={d2:e15}");
		WriteLine($"d1==d2 ? => {d1==d2}");
		WriteLine("\nBy use of the approx function:");
		WriteLine($"d1 approx d2 = approx(d1,d2) = {approx(d1,d2)}");	
	}
}
