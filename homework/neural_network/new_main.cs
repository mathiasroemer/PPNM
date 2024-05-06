using System;
using static System.Math;

public class main{
	public class ann{
		public int n; /* number of hidden neurons */
		public Func<double,double> f = x => x*Exp(-x*x); /* activation function */
		public Func<vector,double> Cost;
		public vector p; /* network parameters */
     		public int steps=0;
		vector x;
		vector y;
		public ann(int n){
			/* constructor */
			this.n = n;
			p = new vector(3*n);
		   	for(int i=0;i<3*n;i++){
				p[i]=1.0;
			}
	   		
			this.Cost = delegate(vector p){
				double sum = 0.0;
				for(int i=0;i<x.size;i++){
					sum += Pow(response(x[i],p)-y[i],2.0);
				}
				
				return sum;
			};
		}
	   public double response(double x, vector p){
	      /* return the response of the network to the input signal x */

		double result = 0;
		for(int i=0;i<3*n;i+=3){
			result += f((x-p[i])/p[i+1])*p[i+2];
		}
		return result;  
	   }
	   public double responseP(double x){
	   	double result = response(x,p);
		return result;
	   }

	   public void train(vector xs,vector ys){
	   /* train the network to interpolate the given table {x,y} */
  		//int iterations = 0;
		/*Func<vector,double> cost_function = delegate(vector p) {
			double result = 0.0;
			for(int i=0;i<x.size;i++){
				result += Pow(response(x[i],p) - y[i],2.0);
			}
			
			double output = result;
		
			iterations++;
			if(iterations % 1000 == 0) {
				Console.Error.WriteLine($"{iterations} {output}");
			}
			return output;
		};*/
		x = xs;
		y = ys;
		p.print("Start guess:");
		p = minimization.newton(Cost,p).Item1;
	   	p.print("Final p:");
	   }
	}
	public static int Main(){
		System.Console.WriteLine("Part A:\n");
		System.Console.WriteLine("First we choose the training function to be Cos(5x-1)*Exp(-x*x)");
		Func<double,double> training_func = x => Cos(5*x-1)*Exp(-x*x);

		int NPoints = 100;
		int NNeurons = 3;
		System.Console.WriteLine($"We then generate {NPoints} random points between -1 and 1 as our input training sample with the points evaluated in the training function as the output");
		vector xs = new vector(NPoints);
		vector ys = new vector(NPoints);

		/*for(int i=0;i<xs.size;i++){
			xs[i] = -1.0+(1.0 - -1.0)*i/(NPoints-1);
			ys[i] = training_func(xs[i]);
		}*/
		for(int i =0;i<NPoints;i++){
        		xs[i] = i*2.0/NPoints-1;
            		ys[i]=training_func(xs[i]);
	        }


		xs.print("X vector [-1,1]:");
		System.Console.WriteLine("");
		ys.print("Y vector training_function(X):");

		System.Console.WriteLine($"We now initiate the neural network with {NNeurons} neurons");
		var nn = new ann(NNeurons);

		System.Console.WriteLine("We then train the network given X and Y");
		nn.train(xs,ys);


		System.Console.Write("Network prediction at 0.2: ");
		System.Console.WriteLine(nn.responseP(0.2));
		System.Console.Write("Real value at 0.2: ");
		System.Console.WriteLine(training_func(0.2));
		
		System.Console.Error.WriteLine("\n");

		for(double i=-1.0;i<1.0+1.0/64.0;i+=1.0/64.0){
			System.Console.Error.WriteLine($"{i} {nn.responseP(i)} {training_func(i)}");
		}

		System.Console.Error.WriteLine("\n");

		for(int i=0;i<xs.size;i++){
			System.Console.Error.WriteLine($"{xs[i]} {ys[i]}");
		}
		return 0;
	}
}
