using System;
using static System.Math;

public class main{
	public class ann{
	   public int n; /* number of hidden neurons */
	   public Func<double,double> f = x => x*Exp(-x*x); /* activation function */
	   public Func<double, double> f_derivative = x => (1.0-2.0*x*x)*Exp(-x*x);
	   public Func<double, double> f_second_derivative = x => 2.0*Exp(-x*x)*x*(2.0*x*x-3.0);
	   public Func<double, double> f_anti_derivative = x => -Exp(-x*x)/2.0;
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
	   public double response_derivative(double x, vector p=null){
	   	if (p==null) p=this.p;

		//F_p(x) = Sum_i f((x-a_i)/b_i) * w_i
		//-> d/dx F_p(x) = Sum_i * f'((x-a_i)/b_i) * (w_i/b_i)
		
		double result = 0;
		for(int i=0;i<p.size;i+=3){
			result += f_derivative((x-p[i])/p[i+1]) * (p[i+2]/p[i+1]);
		}
		return result;
	   }
	   public double response_second_derivative(double x, vector p=null){
	   	if (p==null) p=this.p;

		//d/dx F_p(x) = Sum_i * f'((x-a_i)/b_i) * (w_i/b_i)
		//-> d^2/dx^2 F_p(x) = Sum_i * f''((x-_ai)/b_i) * (w_i/(b_i^2))

		double result = 0;
		for(int i=0;i<p.size;i+=3){
			result += f_second_derivative((x-p[i])/p[i+1]) * (p[i+2]/(p[i+1]*p[i+1]));
		}
		return result;
	   }
	   public double response_anti_derivative(double x, vector p=null){
	   	if (p==null) p=this.p;

		//F_p(x) = Sum_i f((x-a_i)/b_i) * w_i
		//int F_p(x) dx = int Sum_i f((x-a_i)/b_i) * w_i dx
		// = Sum_i w_i * int f((x-a_i)/b_i) dx
		// = Sum_i w_i * int (x-a_i)/b_i)*Exp(-(x-a_i)/b_i)*(x-a_i)/b_i)) dx
		// y=(x-a_i)/b_i, dy/dx = 1/b_i -> dx = dy*b_i
		// -> = Sum_i w_i * b_i * int y*Exp(-y*y) dy
		// = Sum_i w_i * b_i * (-Exp(-y*y)/2)
		// = Sum_i w_i * b_i * (-Exp(-(x-a_i)/b_i * (x-a_i)/b_i)/2)

		double result = 0;
		for(int i=0;i<p.size;i+=3){
			result += p[i+1] * p[i+2] * f_anti_derivative((x-p[i])/p[i+1]);
		}
		return result;

	   }
	   public void train(vector x,vector y){
	   /* train the network to interpolate the given table {x,y} */
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

	   public void train_diff_eq(vector x, vector y,Func<vector,double> phi,double a, double b, vector yc, double alpha=1000, double beta=1000) {
	   /* train the network to interpolate the given table {x,y} */
		Func<vector,double> cost_function = delegate(vector param) {
			
			Func<double,double> phi_integrate = delegate(double k) {
				
				Func<double, vector> diff_vec = delegate(double t) {
					var _temp = new vector(3);
					_temp[0] = this.response_second_derivative(t,param);
					_temp[1] = this.response_derivative(t,param);
					_temp[2] = this.response(t,param);
					return _temp;
				};
				return Math.Pow(phi(diff_vec(k)),2);
			};

			double output = integration.integrate(phi_integrate,a,b).Item1+alpha*Math.Pow(this.response(yc[0],param)-yc[1],2)+beta*Math.Pow(this.response_derivative(yc[2],param)-yc[3],2);
			return output;
		};
		
		(p, steps) = minimization.newton(cost_function,p);
	   }
	}
	public static int Main(){
		System.Console.WriteLine("Part A:\n");
		System.Console.WriteLine("First we choose the training function to be Cos(5x-1)*Exp(-x*x)");
		Func<double,double> training_func = x => Cos(5*x-1)*Exp(-x*x);
		Func<double,double> training_func_derivative = x => Exp(-x*x) * (-5*Sin(5*x - 1) - 2*x*Cos(5*x - 1)); 
		Func<double,double> training_func_second_derivative = x => Exp(-x*x) * ((4*x*x - 27)*Cos(5*x-1)+20*x*Sin(5*x-1)); 

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
			System.Console.Error.WriteLine($"{i} {nn.response(i)} {training_func(i)}");
		}

		System.Console.Error.WriteLine("\n");

		for(int i=0;i<xs.size;i++){
			System.Console.Error.WriteLine($"{xs[i]} {ys[i]}");
		}

		System.Console.Error.WriteLine("\n");

		System.Console.WriteLine("\nPart B:\n");
		System.Console.WriteLine("We then try to calculate the first, second and anti derivative of the training function using the network.");
		
		for(double i=-1.0;i<1.0+1.0/64.0;i+=1.0/64.0){
			System.Console.Error.WriteLine($"{i} {nn.response_derivative(i)} {training_func_derivative(i)}");
		}

		System.Console.Error.WriteLine("\n");

		for(double i=-1.0;i<1.0+1.0/64.0;i+=1.0/64.0){
			System.Console.Error.WriteLine($"{i} {nn.response_second_derivative(i)} {training_func_second_derivative(i)}");
		}

		System.Console.Error.WriteLine("\n");

		for(double i=-1.0;i<1.0+1.0/64.0;i+=1.0/64.0){
			System.Console.Error.WriteLine($"{i} {nn.response_anti_derivative(i)-nn.response_anti_derivative(-1.0)} {integration.integrate(training_func,-1.0,i).Item1}");
		}

		System.Console.WriteLine("\nPart C:\n");

		System.Console.WriteLine("We now want to make the network approximate the solution of a differential equation. For this we test with y´´ = -5y where y(0) = 1 and y´(0) = 1");

		var diffAnn = new ann(NNeurons);
		

		vector yc_initials = new vector(4);
		yc_initials[0] = 0.0;
		yc_initials[1] = 0.0;
		yc_initials[2] = 0.0;
		yc_initials[3] = 1.0;

		Func<vector,double> diff_eq = delegate(vector y){return y[0]+5*y[2];};
		diffAnn.train_diff_eq(null,null,diff_eq,0,2,yc_initials);

		System.Console.Error.WriteLine("\n");

		for(double i=0.0;i<2.0+1.0/64.0;i+=1.0/64.0){
			System.Console.Error.WriteLine($"{i} {diffAnn.response(i)} {Math.Sin(Math.Sqrt(5)*i)/Math.Sqrt(5)}");
		}


		return 0;
	}
}
