class main {
	public static int Main() {
		for(double x=0;x<=5;x+=1.0/8){//we use doubles that can be represented by binary numbers, we can use 5 because it is 101 in binary
			System.Console.WriteLine($"{x} {sfuns.gamma(x)}");
		}

		return 0;
	}
}
