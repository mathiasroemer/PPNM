using System;
using static System.Console;
using static System.Math;

class main{
	public static void Main(string[] args) {
		foreach(var arg in args) {
			var arguments = arg.Split(":");
			if(arguments[0]=="-numbers"){
				var numbers = arguments[1].Split(",");
				foreach (var num in numbers){
					double x = double.Parse(num);
					WriteLine($"{x}: Sin({x}) = {Sin(x)}, Cos({x}) = {Cos(x)}");
				}
			}
		}
	}
}
