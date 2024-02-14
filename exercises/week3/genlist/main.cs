using System;
using static System.Console;

class main{
	public static void Main(){
		var list = new genlist<double[]>();
		char[] delimiters = {' ','\t'};
		var options = StringSplitOptions.RemoveEmptyEntries;
		for(string line = ReadLine();line!=null;line=ReadLine()){
			var numbers = line.Split(delimiters, options);
			int n = numbers.Length;
			var numberList = new double[n];
			for(int i=0;i<n;i++){
				numberList[i] = double.Parse(numbers[i]);
			}
			list.add(numberList);
		}
		for(int i = 0;i<list.size;i++){
			var numbers = list[i];
			foreach(var number in numbers)Write($"{number : 0.00e+00;-0.00e+00} ");
			WriteLine();
		}
		list.remove(0);
		
		WriteLine("Remove element 0:");
		
		for(int i = 0;i<list.size;i++){
			var numbers = list[i];
			foreach(var number in numbers)Write($"{number : 0.00e+00;-0.00e+00} ");
			WriteLine();
		}
		WriteLine("Chain of nodes:");
		list<int> a = new list<int>();
                a.add(1);
                a.add(2);
                a.add(3);
                for( a.start(); a.current != null; a.next()){
                        WriteLine(a.current.item);
                }
	}
}
