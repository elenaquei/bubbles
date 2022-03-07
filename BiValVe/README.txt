26 February 2022

The core library presented here can handle the following cases:
- ODEs:
    - polynomial vector fields
    - computation of single orbit, branch, or 2D-manifold
    - validation of single orbit, branch, or 2D-manifold
    - a posteriori validation of 2D-manifolds
    - automatic blow up for the computation of a Hopf bifurcation
    - saddle node validation w.r.t. any variable
    - Hopf bifurcation validation (with blow up approach)
 following the papers: 
        1) van den Berg JB, Queirolo E. A general framework for validated continuation of periodic orbits in systems of polynomial ODEs. Journal of Computational Dynamics. 2021;8(1):59.
        2) van den Berg JB, Lessard JP, Queirolo E. Rigorous verification of Hopf bifurcations via desingularization and continuation. SIAM Journal on Applied Dynamical Systems. 2021;20(2):573-607.
        3) Church KM, Queirolo E. Computer-assisted proofs of Hopf bubbles and degenerate Hopf bifurcations, Submitted Feb 2022

- DDEs:
    - polynomial vector fields with non-polynomial scalar equations
    - computation of 2D-manifolds 
    - validation of 2D-manifold
    - a posteriori validation of 2D-manifolds
following:
        Church KM, Queirolo E. Computer-assisted proofs of Hopf bubbles and degenerate Hopf bifurcations, Submitted Feb 2022

Initialise the library by running
    start_BiValVe
Include all examples by running
    load_examples
Do not move the working directory to any folder.

For an introduction to the library, the folder "branch ODE templates" present some easy examples.
More examples are presented in "branch Hopf examples".
All examples for 
        Church KM, Queirolo E. Computer-assisted proofs of Hopf bubbles and degenerate Hopf bifurcations, Submitted Feb 2022
can be found in the folders "simplex ODE Hopf examples" and "simplex DDE Hopf examples", or can be run with "figure_generation".
Warning: 2D-manifold examples and DDE examples require significant memory.

