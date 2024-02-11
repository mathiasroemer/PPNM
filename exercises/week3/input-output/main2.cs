using System;
using static System.Console;
using static System.Math;

class main2 {
	public static void Main(){
		char[] split_delimiters = {' ', '\t','\n'};
		var split_options = StringSplitOptions.RemoveEmptyEntries;
		for(string line = ReadLine(); line!=null; line=ReadLine()){
			var numbers = line.Split(split_delimiters,split_options);
			foreach(var num in numbers) {
				double x = double.Parse(num);
				Error.WriteLine($"{x}: Sin({x}) = {Sin(x)}, Cos({x}) = {Cos(x)}");
			}
		}
	}	
}
