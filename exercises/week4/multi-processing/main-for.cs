class main {
	
	public static int Main(string[] args){
		int nterms = (int)1e8;
		
		System.Console.Write("The sum variable is not local giving a race condition making the threads trying to access a shared variable which will take more time. Additionally will there be some inaccuracy to the final result since different double are added to eachother.\n");
	
		foreach(string arg in args){
			var words = arg.Split(":");
			if(words[0]=="-nterms") nterms=(int)double.Parse(words[1]);
		}	
		System.Console.Write($"Main: nterms={nterms}\n");


		double sum=0;
		System.Threading.Tasks.Parallel.For(1,nterms+1,delegate(int i){sum+=1.0/i;});

		System.Console.Write($"Main: Total sum = {sum}\n");
		return 0;
	}//Main
}//main
