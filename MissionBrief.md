# What is the point of this project.

The aim is to simulate various random non-commutative geometries. Specifically, to simulate any finite spectral triple as defined in the following paper: [Matrix geometries and fuzzy spaces as finite spectral triples - John W. Barrett](https://arxiv.org/abs/1502.05383).

Why would we want to do this? This is explained in the paper linked above in more detail, or you can read the introduction to my [PhD thesis](https://www.pauldruce.co.uk/assets/PDFs/Paul%20Druce's%20PhD%20Thesis%202020.pdf), but here is a brief over view from a physics perspective.
Noncommutative geometry is an attempt to replace the ordinary geometry that you may be familiar with, with a new framework that places better with quantum mechanics.
Image you want to include the notion of an uncertainty principle into your geometry.
So instead of the usual Heisenberg uncertainty principle, which states you can't know the position and momentum of a quantum object in a specific direction (i.e. x-coordinate and x-momentum) with perfect accuracy.
The noncommutative geometry approach is to introduct the idea that you can't know all the coordinates simultaneusly to perfect precision.
This is the starting idea of noncommutative geometry, and a number of different schemes have been developed to enact this desired uncertainity.

One of the approached is that we change the space of continuous functions over a space to a new space that no longer is commutative.
For instance, the x-position and the y-position are functions of a space, i.e. $\hat{x}(p)= p_x$, where $p$ is a point on the space.
So if we now say that $[\hat{x}, \hat{y}] \neq 0$, then this introduces the desired uncertainty principle between the $x$ and $y$ coordinates.

So the scheme that the paper and thesis above are based upon, is the idea that the space of functions of a geometry are replaced with a finite matrix algebra.
This has deeper justifications (again see my thesis or the paper), and there exists a complete characterisation of the allowed noncommutative geometries with such a matrix algebra for the function space.

It is over these characterisations that the simulation is performed over.

The aim is to find non-commutative geometries that appear to have similar properties as ordinary manifold.
