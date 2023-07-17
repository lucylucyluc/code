In the maincode, there are five code for the paper "Efficient estimation of Cox Model with Random Change Point".

In the first code, "set1_exp",  we assume that the random change point is from a exponential distribution under the Cox Model with Covariate-free Random Change Point  without covariate Z. 
while the random change point is assumed to follow a normal distribution with mean mu and variance sigma^2 in code "set_norm". And the form of cox model is the same as set1_exp.

The code "set2" is for setting 2 in this paper, where the random change point also form a normal distribution. The difference between this code and set1_norm is that the form of the cox model contain covariate Z.

The code "set3" is for setting 3 in this paper, where the random change point has a linear structure with an error. And the error form a normal distribution. In the cox model, covariate Z = ∅.

The last code, "realdata2" is for the breast dataset under Cox Model with Covariate-dependent Random Change Point.
Note that in the real data analysis, the calculation process of Cox Model with Covariate-free Random Change Point can be implement by code set2.

In all code, the "dat" is data, which is a list. The "par0" is the initial value that needs to be input, "k" is the number of inner knots of spline, and "intebound" include the lower and upper bounds of integral.