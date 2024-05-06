using System;
using static System.Math;

public class main{
	public class ann{
	   public int n; /* number of hidden neurons */
	   public Func<double,double> f = x => x*Exp(-x*x); /* activation function */
	   public vector p; /* network parameters */
	   public int steps;
	   public ann(int n){
		/* constructor */
		this.n = n;
		p = new vector(3*n);
	   	for (int i=0;i<p.size;i++){
			p[i]=1.0;
		}
	   }
	   public double response(double x,vector p=null){
	      /* return the response of the network to the input signal x */
		if (p==null) p=this.p;

		double result = 0;
		for(int i=0;i<p.size;i+=3){
			result += f((x-p[i])/p[i+1])*p[i+2];
		}
		return result;  
	   }
	   public void train(vector x,vector y){
	   /* train the network to interpolate the given table {x,y} */
  		int iterations = 0;
		Func<vector,double> cost_function = delegate(vector param) {
			double result = 0.0;
			for(int i=0;i<x.size;i++){
				result += Pow(response(x[i],param) - y[i],2.0);
			}
			
			double output = result/x.size;
		
			return output;
		};
		
		(p, steps) = minimization.newton(cost_function,p);
	   }
	}
	public static int Main(){
		System.Console.WriteLine("Part A:\n");
		System.Console.WriteLine("First we choose the training function to be Cos(5x-1)*Exp(-x*x)");
		Func<double,double> training_func = x => Math.Cos(5*x-1)*Math.Exp(-x*x);
		//Func<double,double> training_func = x => Pow(x,2);

		int NPoints = 100;
		int NNeurons = 3;
		System.Console.WriteLine($"We then generate {NPoints} random points between -1 and 1 as our input training sample (X) with the points evaluated in the training function as the output (Y)");
		vector xs = new vector(NPoints);
		vector ys = new vector(NPoints);

		for(int i=0;i<xs.size;i++){
			xs[i] = -1.0+(1.0 - -1.0)*i/(NPoints-1);
			ys[i] = training_func(xs[i]);
		}

		System.Console.WriteLine($"We then initiate the neural network with {NNeurons} neurons");
		var nn = new ann(NNeurons);

		System.Console.WriteLine("The initial network parameters are:");
		nn.p.print();
		
		System.Console.WriteLine("\nWe then train the network given X and Y\n");
		nn.train(xs,ys);

		System.Console.WriteLine("The network parameters then becomes after training:");
		nn.p.print();

		System.Console.Write("\n\nNetwork prediction at 0.2: ");
		System.Console.WriteLine(nn.response(0.2));
		System.Console.Write("Real value at 0.2: ");
		System.Console.WriteLine(training_func(0.2));
		
		System.Console.Error.WriteLine("\n");

		for(double i=-1.0;i<1.0+1.0/64.0;i+=1.0/64.0){
			System.Console.Error.WriteLine($"{i} {nn.response(i)}");
		}

		System.Console.Error.WriteLine("\n");

		for(int i=0;i<xs.size;i++){
			System.Console.Error.WriteLine($"{xs[i]} {ys[i]}");
		}
		return 0;
	}
}
