using System.Linq;

class main {
	
	public static int Main(string[] args){
		int nterms = (int)1e8;
	
		System.Console.Write($"The result is now more correct but is still slower, once again because of the race condition.\n");

		foreach(string arg in args){
			var words = arg.Split(":");
			if(words[0]=="-nterms") nterms=(int)double.Parse(words[1]);
		}	
		System.Console.Write($"Main: nterms={nterms}\n");

		var sum = new System.Threading.ThreadLocal<double>( ()=>0, trackAllValues:true);
		System.Threading.Tasks.Parallel.For( 1, nterms+1, (int i)=>sum.Value+=1.0/i );
		double totalsum=sum.Values.Sum();

		System.Console.Write($"Main: Total sum = {totalsum}\n");
		return 0;
	}//Main
}//main
