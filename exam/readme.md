Cubic (sub-)spline for data with derivatives has been implemented.
This has been done by solving for the three coefficients b_i, c_i and d_i. They were solved by substituting the spline equation into the three condition equations. The three coefficients were then found to be almost exactly equal to the coefficients given for the Akima subspline in the book. The only difference was that the spline derivatives should be replaced with the data derivatives.

In the svg file "Out.SignSpline.svg" can a plot similar to the one provided in the book be seen. A sign function has been set up where a couple of points has then been given to the subspline together with their derivatives. The resulting subspline interpolation looks much like the one in the book. For comparison reasons, the cubic spline found in the Spline homework has also been utilized.

Going further to the svg file "Out.CosDerivAndInteg.svg", the subspline derivative and integral has been used together with the subspline. This is done on a Cosine function, as the resulting derivative and integral make it easy to conclude that the derivative and integration methods works.

In the "Out.SecondDerivative.svg" file I attempted to make the second derivative continuous by increasing the order of the spline from 3 to 4. The coefficients b_i, c_i and d_i were found to be exactly equal to the ones found in the order 3 case when substituting the spline into the conditions. In order to get a continuous second derivative I tried using the condition "S_i''(x_i+1) = S_i+1''(x_i+1)" to find the e_i coefficients.
The resulting e_i coefficients was found such that e_i+1 required e_i. In order to then calculate the e_i's I tried to average the resulting forward and backward recursion coefficients. As seen in the plot, this did not work as it was supposed to.

In general I do believe that the original order 3 subspline works as intended and do show an improved job on the sign function compared to the cubic spline interpolation. I also believe that the first derivative and the anti derivative methods works as intended. I am a little disappointed that I could not make the second derivative continuous but everything else seems to function.

Self-evaluation: (9/10)
