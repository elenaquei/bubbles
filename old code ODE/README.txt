3 June 2020

This program validates a branch of periodic orbits of a given ODE system.
A first solution is computed, then the branch is continued by one step and the segment 
in-between is validated. If the validation succeeds, a new step is performed, if the
validation failed, the segment is recomputed doing a smaller step size (and eventually
adding nodes to the solution).
The latest extension to the core library allows the user to validate both Hopf bifurcations and saddle node bifurcations.

This library consists of a variety of functions whose aim is to validate a point wise solution, or a segment or a branch of solutions, considering just periodic solutions to polynomial ODEs. The added functionalities use the results of the validation to prove the existence (or absence) of saddle nodes or Hopf bifurcations.

The more important functions you should be concerned with are

continuation (a simple function to help you start your continuation proof, require a non-square problem, including scalar equations, and a single solution as in instance of Xi_vector) - it has many other optional parameters, check its help for further information

radii_polynomails_cont_new (a more involved function validating a segment between two approximate solutions)

continuation_Hopf (similar to continuation, uses some information on the numerical Hopf bifurcation to set up the validation of the branch of periodic orbits)

algebraic_hopf (validate the existence of a Hopf bifurcation considering only the eigenvalue problem of the equilibrium at the Hopf bifurcation. The outputs of this function are useful to start continuation_Hopf)

validation_saddle_or_not (checks if a segment of the branch is undergoing a saddle node bifurcation and validates the existence of absence of the saddle node)

In order to help you approach this two main functions, a couple of helpers functions are there for you:

from_string_to_polynomial_coef is a function that takes as input a string defining a vector field and construct an instance of the class polynomial_coef that corresponds to that string (please view the help function to write the string according to the specified format)

time_series2Xi_vec is a function that takes as input a period and a time series and returns an instance of the class Xi_vector to work with

fancy_scalar_condition takes as input the constructed Xi_vector and creates an instance of the class scalar_eq corresponding to fixing the the shift of your Xi_vector (i.e. <x',x>= 0)

Remark: this specific code is built in such a way to need Intlab, therefore Intlab must be in the path. Furthermore, if Intlab hasn't been initialised yet, the code will initialise it, therefore giving an error very soon on if Intlab is not running.

For a first approach to the library, check out the template folder, including the files template_vdP, a script showing you how to use the functions just highlighted to get a point wise validation, template_lor_cont shows you how to validate a branch of solution in the simplest way.
A last example provided is template_vdP_cont, that shows the difference for running continuation in the van der Pol system.


If you are interested in generating the figures presented in the second accompanying paper "Rigorous verification of Hopf bifurcations via desingularization and continuation" by van den Berg, Lessard and Queirolo, the folder "figures_articleHopf" has all the codes specifically pertaining to the article and "all_figure" creates all the figures in order. The code runs in some 6 hours total.