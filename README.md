# geod
geodesic integration for chaos searches


This is a collection of code made to perform numerical integration and then do analysis of that data.

The integrator is written in C and can be found in geod_src.
It is designed to integrate trajectories in arbitrary 6-dimensional phase spaces, and is highly modular,
so that adding a new equation of motion just means adding one more C source file.

It was written at UIUC in order to do very high-resolution integrations of geodesics to look for chaos
in modified theories of gravity.  To that end, it is an implementation of an RKF7(8) integrator, meaning
it is accurate to O(h^8) for timestep h.  This means that you can take relatively large timesteps and still
have very high precision.

The analysis software is all written in Python.  This deals with constructing Poincare maps, calculating
rotation numbers, and creating helpful functions for dealing with large directories of big data files.
There's also some code designed to make it easy to run large grids of initial conditions on a cluster.
In addition, there's a Python front-end to call the C integrator, which is really nice to work with in
a notebook.

At the moment, the code is rather messy (the C in particular is filled with unused functions and debugging
print statements), and so I'm not ready to put documentation up yet.  However, it is a fast and good integrator
with wide applications, and a nice command-line interface.  If you want to use it, let me know and I'll walk you
through it.
