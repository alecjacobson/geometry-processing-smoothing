# Geometry Processing - Smoothing

> **To get started:** Clone this repository then issue
> 
>     git clone --recursive http://github.com/[username]/geometry-processing-smoothing.git
>

## Installation, Layout, and Compilation

See
[introduction](http://github.com/alecjacobson/geometry-processing-introduction).

## Execution

Once built, you can execute the assignment from inside the `build/` by running
on a given mesh with given scalar field (in
[dmat](http://libigl.github.io/libigl/file-formats/dmat.html) format).

    ./smoothing [path to mesh.obj] [path to data.dmat]

or to load a mesh with phony noisy data use:

    ./smoothing [path to mesh.obj] n

or to load a mesh with smooth z-values as data (for mesh smoothing only):

    ./smoothing [path to mesh.obj]


## Background

In this assignment we will explore how to smooth a data _signal_ defined over a
curved surface. The data _signal_ may be a scalar field defined on a static
surface: for example, noisy temperatures on the surface of an airplane.
Smoothing in this context can be understood as [data
denoising](https://en.wikipedia.org/wiki/Noise_reduction):

![](images/plane-smooth-data.gif)

The _signal_ could also be the geometric positions of the surface itself. In
this context, smoothing acts also affects the underlying geometry of the
domain. We can understand this operation as [surface
fairing](https://en.wikipedia.org/wiki/Surface_fairng):

![](images/plane-smooth-geometry.gif)

### Flow-based formulation

In both cases, the underlying mathematics of both operations will be very
similar. If we think of the signal as undergoing a _flow_ toward a smooth
solution over some phony notion of "time", then the governing partial
differential equation we will start with sets the change in signal value <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/>
over time [proportional
to](https://en.wikipedia.org/wiki/Proportionality_(mathematics)) the
[Laplacian](https://en.wikipedia.org/wiki/Laplace_operator) of the signal <img src="./tex/520957eebcbc4341addf55ec1a6ef0ab.svg?invert_in_darkmode" align=middle width=23.10894465pt height=22.4657235pt/>
(for now, roughly the second derivative of the signal as we move _on_ the
surface):

<p align="center"><img src="./tex/eb0a14f1824b0e6901ad51c2b68e8ee7.svg?invert_in_darkmode" align=middle width=80.20514865pt height=33.8120871pt/></p>


where the scalar parameter <img src="./tex/1b109d8b4484cf614f27126d788c510e.svg?invert_in_darkmode" align=middle width=9.58908225pt height=22.8310566pt/> controls the rate of smoothing.

When the signal is the surface geometry, we call this a [geometric
flow](https://en.wikipedia.org/wiki/Geometric_flow).

There are various ways to motivate this choice of flow for
data-/geometry-smoothing. Let us consider one way that will introduce the
Laplacian as a form of local averaging.

Given a noisy signal <img src="./tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode" align=middle width=9.81741585pt height=22.8310566pt/>, intuitively we can _smooth_ <img src="./tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode" align=middle width=9.81741585pt height=22.8310566pt/> by averaging every
value with its neighbors' values. In continuous math, we might write that the
smoothed value <img src="./tex/193e99dd587ecc23cdd8853e4cd46a76.svg?invert_in_darkmode" align=middle width=32.172822pt height=24.657534pt/> at any point on our surface <img src="./tex/191fb656c7f4f20997cf240f7af26caf.svg?invert_in_darkmode" align=middle width=40.5704805pt height=22.5570873pt/> should be equal to
the average value of some small
[ball](https://en.wikipedia.org/wiki/Ball_(mathematics)) of nearby points:

<p align="center"><img src="./tex/a69ec5fa9866f2be8caf1efd6fdc8d11.svg?invert_in_darkmode" align=middle width=210.76788645pt height=40.2286731pt/></p>


If the ball <img src="./tex/06e8e88b4ac31efcc03ecf253f9609ba.svg?invert_in_darkmode" align=middle width=36.05595345pt height=24.657534pt/> is small, then we will have to repeat this averaging many
times to see a global smoothing effect. Hence, we can write that the current
value <img src="./tex/2781fb506c10958a89548aa070f15594.svg?invert_in_darkmode" align=middle width=14.37606555pt height=26.0859621pt/> _flows_ toward smooth solution by small steps <img src="./tex/c6bcca8f259a226ca0748027158124bd.svg?invert_in_darkmode" align=middle width=13.86416955pt height=22.8310566pt/> in time:

<p align="center"><img src="./tex/e8b064f3dd720aa90964378e9d8b8f47.svg?invert_in_darkmode" align=middle width=243.3853983pt height=40.2286731pt/></p>


Subtracting the current value <img src="./tex/4813491b0a749868771cab60ce88649f.svg?invert_in_darkmode" align=middle width=37.9605105pt height=26.0859621pt/> from both sides and introducing a
flow-speed parameter <img src="./tex/1b109d8b4484cf614f27126d788c510e.svg?invert_in_darkmode" align=middle width=9.58908225pt height=22.8310566pt/> we have a flow equation
describing the change in value as an integral of relative values:

<p align="center"><img src="./tex/32bab51d0ad9d94ce1fe47d67e452dbb.svg?invert_in_darkmode" align=middle width=271.1100612pt height=40.3902642pt/></p>


For harmonic functions, <img src="./tex/bb23bb09b67bb167730c76ff706043de.svg?invert_in_darkmode" align=middle width=53.24578545pt height=22.4657235pt/>, this integral becomes zero in the limit as the
radius of the ball shrinks to zero via satisfaction of the
[mean value
theorem](https://en.wikipedia.org/wiki/Harmonic_function#The_mean_value_property).
It follows for a non-harmonic <img src="./tex/74a93b9b98af7c7b321f354c4fcf8ea2.svg?invert_in_darkmode" align=middle width=53.24578545pt height=22.8310566pt/> this integral is equal to the Laplacian
of the <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/>, so we have arrived at our flow equation:

<p align="center"><img src="./tex/1837a3ba67cc7b4f946656323866bfcb.svg?invert_in_darkmode" align=middle width=384.5539698pt height=40.3902642pt/></p>


### Energy-based formulation

Alternatively, we can think of a single smoothing operation as the solution to
an energy minimization problem. If <img src="./tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode" align=middle width=9.81741585pt height=22.8310566pt/> is our noisy signal over the surface,
then we want to find a signal <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/> such that it simultaneously minimizes its
difference with <img src="./tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode" align=middle width=9.81741585pt height=22.8310566pt/> and minimizes its variation over the surface:

<p align="center"><img src="./tex/2103cd36ae87c74dbd5cf1297f2bb1da.svg?invert_in_darkmode" align=middle width=419.77996005pt height=49.05031725pt/></p>


where again the scalar parameter <img src="./tex/1b109d8b4484cf614f27126d788c510e.svg?invert_in_darkmode" align=middle width=9.58908225pt height=22.8310566pt/> controls the rate of smoothing. This
energy-based formulation is equivalent to the flow-based formulation.
Minimizing these energies is identical to stepping forward one temporal unit in
the flow.

#### Calculus of variations
In the smooth setting, our energy <img src="./tex/84df98c65d88c6adf15d4645ffa25e47.svg?invert_in_darkmode" align=middle width=13.0821966pt height=22.4657235pt/> is a function that measures scalar value
of a given function _u_, making it a
[functional](https://en.wikipedia.org/wiki/Functional_(mathematics)). To
understand how to _minimize_ a functional with respect to an unknown function,
we will need concepts from the [calculus of
variations](https://en.wikipedia.org/wiki/Calculus_of_variations).

We are used to working with minimizing quadratic _functions_ with respect to a
discrete set of variables, where the minimum is obtained when the gradient of
the energy with respect to the variables is zero.

In our case, the functional <img src="./tex/36f9bbd82b66589fd9772b7889d95932.svg?invert_in_darkmode" align=middle width=35.27788275pt height=24.657534pt/> is quadratic in <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/> (recall that the
[gradient operator](https://en.wikipedia.org/wiki/Gradient) <img src="./tex/289bc3c68d369342d98da6ff43ec5b9f.svg?invert_in_darkmode" align=middle width=13.69867125pt height=22.4657235pt/> is a linear
operator). The function <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/> that minimizes <img src="./tex/36f9bbd82b66589fd9772b7889d95932.svg?invert_in_darkmode" align=middle width=35.27788275pt height=24.657534pt/> will be obtained when any
small change or _variation_ in <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/> has no change on the energy values. To
create a small change in a function <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/> we will add another function <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.5578603pt height=14.1552444pt/> times
a infinitesimal scalar <img src="./tex/7a676caa07476cecc1fef337bd6a860e.svg?invert_in_darkmode" align=middle width=6.6723921pt height=14.1552444pt/>. If <img src="./tex/36f9bbd82b66589fd9772b7889d95932.svg?invert_in_darkmode" align=middle width=35.27788275pt height=24.657534pt/> is minimized for a function <img src="./tex/31fae8b8b78ebe01cbfbe2fe53832624.svg?invert_in_darkmode" align=middle width=12.21084645pt height=14.1552444pt/> and we
are given another arbitrary function <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.5578603pt height=14.1552444pt/>, then let us define a function new
function 

<p align="center"><img src="./tex/1676aad3b697e1ce24e6fdf9e3b80555.svg?invert_in_darkmode" align=middle width=448.4368152pt height=37.3519608pt/></p>

where we observe that <img src="./tex/eb62325152b7dd734820edcdd33d2950.svg?invert_in_darkmode" align=middle width=9.79454355pt height=22.8310566pt/> is quadratic in <img src="./tex/7a676caa07476cecc1fef337bd6a860e.svg?invert_in_darkmode" align=middle width=6.6723921pt height=14.1552444pt/>.

Since <img src="./tex/21d2cb0fed44339d6eb9fd1c6788522f.svg?invert_in_darkmode" align=middle width=38.0784558pt height=24.657534pt/> is minimal then <img src="./tex/eb62325152b7dd734820edcdd33d2950.svg?invert_in_darkmode" align=middle width=9.79454355pt height=22.8310566pt/> is minimized when <img src="./tex/7a676caa07476cecc1fef337bd6a860e.svg?invert_in_darkmode" align=middle width=6.6723921pt height=14.1552444pt/> is zero, and if <img src="./tex/eb62325152b7dd734820edcdd33d2950.svg?invert_in_darkmode" align=middle width=9.79454355pt height=22.8310566pt/> is
minimal at <img src="./tex/00e14c64b601ac983b4364aa60fa1dc6.svg?invert_in_darkmode" align=middle width=36.80923125pt height=21.1872144pt/>, then the derivative of <img src="./tex/eb62325152b7dd734820edcdd33d2950.svg?invert_in_darkmode" align=middle width=9.79454355pt height=22.8310566pt/> with respect <img src="./tex/7a676caa07476cecc1fef337bd6a860e.svg?invert_in_darkmode" align=middle width=6.6723921pt height=14.1552444pt/> must be zero:


<p align="center"><img src="./tex/f09d1063c8d652ad8acdd78ea9db81f3.svg?invert_in_darkmode" align=middle width=695.17788615pt height=224.75072565pt/></p>


The choice of "test" function <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.5578603pt height=14.1552444pt/> was arbitrary, so this must hold for any
(reasonable) choice of <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.5578603pt height=14.1552444pt/>:

<p align="center"><img src="./tex/f918be183bda0a2a6d22fe2690f24395.svg?invert_in_darkmode" align=middle width=278.4682362pt height=37.3519608pt/></p>


It is difficult to claim much about <img src="./tex/31fae8b8b78ebe01cbfbe2fe53832624.svg?invert_in_darkmode" align=middle width=12.21084645pt height=14.1552444pt/> from this equation directly because
derivatives of <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.5578603pt height=14.1552444pt/> are still involved. We can _move_ a derivative from <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.5578603pt height=14.1552444pt/> to a
<img src="./tex/31fae8b8b78ebe01cbfbe2fe53832624.svg?invert_in_darkmode" align=middle width=12.21084645pt height=14.1552444pt/> by applying [Green's first
identity](https://en.wikipedia.org/wiki/Green's_identities):

<p align="center"><img src="./tex/04c2c8b27440c9aea318748e3a1a9ceb.svg?invert_in_darkmode" align=middle width=402.7609938pt height=37.3519608pt/></p>

where we choose to _ignore_ the boundary terms (for now) or equivalently we
agree to work on _closed_ surfaces <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.5022269pt height=22.5570873pt/>.

Since this equality must hold of _any_ <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.5578603pt height=14.1552444pt/> let us consider functions that are
little ["blips"](https://en.wikipedia.org/wiki/Bump_function) centered at any
arbitrary point <img src="./tex/191fb656c7f4f20997cf240f7af26caf.svg?invert_in_darkmode" align=middle width=40.5704805pt height=22.5570873pt/>. A function <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.5578603pt height=14.1552444pt/> that is one at <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711605pt height=14.6118786pt/> and quickly
decays to zero everywhere else. To satisfy the equation above at <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711605pt height=14.6118786pt/> with this
blip <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.5578603pt height=14.1552444pt/> we must have that:

<p align="center"><img src="./tex/e94a6a9790d776ead8b76d74b10cc90f.svg?invert_in_darkmode" align=middle width=172.38955305pt height=16.438356pt/></p>

The choice of <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711605pt height=14.6118786pt/> was arbitrary so this must hold _everywhere_.

Because we invoke _variations_ to arrive at this equation, we call the
_energy-based_ formulation a _variational formulation_.

### Implicit smoothing iteration

Now we understand that the flow-based formulation and the variational
formulation lead to the same system, let us concretely write out the implicit
smoothing step. 

Letting <img src="./tex/eb367462f49f6165b7702489cd3ad4ea.svg?invert_in_darkmode" align=middle width=48.51977955pt height=26.7617526pt/> we compute a new smoothed function <img src="./tex/74341d0452aa58e0f8653b6058b1f3f6.svg?invert_in_darkmode" align=middle width=31.01998515pt height=26.7617526pt/> given the
current solution <img src="./tex/2781fb506c10958a89548aa070f15594.svg?invert_in_darkmode" align=middle width=14.37606555pt height=26.0859621pt/> by solving the _linear_ system of equations:

<p align="center"><img src="./tex/ca2a12218ffab9d11bc2166ed9ad9083.svg?invert_in_darkmode" align=middle width=257.79280725pt height=18.3123831pt/></p>

where <img src="./tex/2b41e5149bbd9bd2a0b7170924be9b80.svg?invert_in_darkmode" align=middle width=13.69867125pt height=22.8310566pt/> is the [identity
operator](https://en.wikipedia.org/wiki/Identity_function). In the discrete
case, we will need discrete approximations of the <img src="./tex/2b41e5149bbd9bd2a0b7170924be9b80.svg?invert_in_darkmode" align=middle width=13.69867125pt height=22.8310566pt/> and <img src="./tex/545994479b8eb0934d80a97627e63f82.svg?invert_in_darkmode" align=middle width=13.69867125pt height=22.4657235pt/>
operators.

## Discrete Laplacian

There are many ways to derive a discrete
approximation of the Laplacian <img src="./tex/545994479b8eb0934d80a97627e63f82.svg?invert_in_darkmode" align=middle width=13.69867125pt height=22.4657235pt/> operator on a triangle mesh using:

 - [finite volume method](https://en.wikipedia.org/wiki/Finite_volume_method),
    - "The solution of partial differential equations by means of electrical
      networks"  [MacNeal 1949, pp. 68],
    - "Discrete differential-geometry operators for triangulated
      2-manifolds"  [Meyer et al. 2002],
    - _Polygon mesh processing_ [Botsch et al. 2010],
 - [finite element
   method](https://en.wikipedia.org/wiki/Finite_element_method),
    - "Variational methods for the solution of problems of equilibrium and
      vibrations" [Courant 1943],
    - _Algorithms and Interfaces for Real-Time Deformation of 2D and 3D Shapes_
      [Jacobson 2013, pp. 9]
 - [discrete exterior
   calculus](https://en.wikipedia.org/wiki/Discrete_exterior_calculus)
    - _Discrete Exterior Calculus_ [Hirani 2003, pp. 69]
    - _Discrete Differential Geometry: An Applied Introduction_ [Crane 2013,
      pp. 71]
 - [gradient of surface area \Rightarrow  mean curvature
   flow](https://en.wikipedia.org/wiki/Mean_curvature_flow)
    - "Computing Discrete Minimal Surfaces and Their Conjugates" [Pinkall &
      Polthier 1993]

All of these techniques will produce the _same_ sparse _Laplacian matrix_ <img src="./tex/741ed3ef6a246c4afd8317583b95f093.svg?invert_in_darkmode" align=middle width=69.85918335pt height=26.1773094pt/>
for a mesh with <img src="./tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode" align=middle width=9.86687625pt height=14.1552444pt/> vertices. 

### Finite element derivation of the discrete Laplacian

We want to approximate the Laplacian of a function <img src="./tex/520957eebcbc4341addf55ec1a6ef0ab.svg?invert_in_darkmode" align=middle width=23.10894465pt height=22.4657235pt/>. Let us consider <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/> to
be [piecewise-linear](https://en.wikipedia.org/wiki/Piecewise_linear_function)
represented by scalar values at each vertex, collected in <img src="./tex/1a693d20e239e8c342dd2a059e9770c6.svg?invert_in_darkmode" align=middle width=50.5915674pt height=22.6483917pt/>.

Any piecewise-linear function can be expressed as a sum of values at mesh
vertices times corresponding piecewise-linear basis functions  (a.k.a hat
functions, <img src="./tex/965b9779f7d35773b987d1080ab428e9.svg?invert_in_darkmode" align=middle width=15.40433235pt height=14.1552444pt/>):


<p align="center"><img src="./tex/e0b10549ef06a4ad5f8b1880b7d31832.svg?invert_in_darkmode" align=middle width=330.96304065pt height=105.7206117pt/></p>



<p align="center"><img src="./tex/0ae2ff8f9458d7ebf324995914023135.svg?invert_in_darkmode" align=middle width=492.09227925pt height=69.0417981pt/></p>


![](images/hat-function.png)

By plugging this definition into our smoothness energy above, we have discrete
energy that is quadratic in the values at each mesh vertex:


<p align="center"><img src="./tex/f708fca18bc40e08e6d26508a0f46968.svg?invert_in_darkmode" align=middle width=465.1655811pt height=169.13386875pt/></p>



By defining <img src="./tex/965b9779f7d35773b987d1080ab428e9.svg?invert_in_darkmode" align=middle width=15.40433235pt height=14.1552444pt/> as piecewise-linear hat functions, the values in the system
matrix <img src="./tex/58e56656e6e782a9941daf494e8471eb.svg?invert_in_darkmode" align=middle width=21.9426504pt height=22.4657235pt/> are uniquely determined by the geometry of the underlying mesh.
These values are famously known as _cotangent weights_. "Cotangent"
because, as we will shortly see, of their trigonometric formulae and "weights"
because as a matrix <img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979855pt height=22.5570873pt/> they define a weighted [graph
Laplacian](https://en.wikipedia.org/wiki/Laplacian_matrix) for the given mesh.
Graph Laplacians are employed often in geometry processing, and often in
discrete contexts ostensibly disconnected from FEM. The choice or manipulation
of Laplacian weights and subsequent use as a discrete Laplace operator has been
a point of controversy in geometry processing research (see "Discrete laplace
operators: no free lunch" [Wardetzky et al. 2007]).


We first notice that <img src="./tex/7d9f97a055843950ec03bc359085c06e.svg?invert_in_darkmode" align=middle width=29.1030036pt height=22.4657235pt/> are constant on each triangle, and only nonzero on
triangles incident on node <img src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.6632257pt height=21.6830097pt/>. For such a triangle, <img src="./tex/04fdc9ebe2a4f15ffdd732d0abd189e1.svg?invert_in_darkmode" align=middle width=18.1520988pt height=22.4657235pt/>, this <img src="./tex/7d9f97a055843950ec03bc359085c06e.svg?invert_in_darkmode" align=middle width=29.1030036pt height=22.4657235pt/> points
perpendicularly from the opposite edge <img src="./tex/b95c2b0aab2482e5bebd25332a4bbde0.svg?invert_in_darkmode" align=middle width=12.3050367pt height=14.1552444pt/> with inverse magnitude equal to
the height <img src="./tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode" align=middle width=9.4711155pt height=22.8310566pt/> of the triangle treating that opposite edge as base:

<p align="center"><img src="./tex/4148b1898d976255af11ef4509e668c1.svg?invert_in_darkmode" align=middle width=142.70184225pt height=34.7253258pt/></p>

where <img src="./tex/2dfa893f7bf3c5c01ceb4788da782da0.svg?invert_in_darkmode" align=middle width=13.3152525pt height=14.6118786pt/> is the edge <img src="./tex/b95c2b0aab2482e5bebd25332a4bbde0.svg?invert_in_darkmode" align=middle width=12.3050367pt height=14.1552444pt/> as a vector and <img src="./tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode" align=middle width=12.32879835pt height=22.4657235pt/> is the area of the triangle.

![Left: the gradient <img src="./tex/461fd2586080bcfe7ac91515a525affe.svg?invert_in_darkmode" align=middle width=29.1030036pt height=22.4657235pt/> of a hat function <img src="./tex/965b9779f7d35773b987d1080ab428e9.svg?invert_in_darkmode" align=middle width=15.40433235pt height=14.1552444pt/> is piecewise-constant and
points perpendicular to opposite edges. Right: hat function gradients <img src="./tex/461fd2586080bcfe7ac91515a525affe.svg?invert_in_darkmode" align=middle width=29.1030036pt height=22.4657235pt/>
and <img src="./tex/e7299099834dc07da267043c72b8955e.svg?invert_in_darkmode" align=middle width=30.556614pt height=22.4657235pt/> of neighboring nodes meet at angle <img src="./tex/af94b55c890d07e5e8a5cbccff7647c0.svg?invert_in_darkmode" align=middle width=81.4134915pt height=22.8310566pt/>.](images/hat-function-gradient.png)

Now, consider two neighboring nodes <img src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.6632257pt height=21.6830097pt/> and <img src="./tex/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710417pt height=21.6830097pt/>, connected by some edge
<img src="./tex/902505b1f10837bf329ffbeb68762c74.svg?invert_in_darkmode" align=middle width=19.41976245pt height=14.6118786pt/>. Then <img src="./tex/7d9f97a055843950ec03bc359085c06e.svg?invert_in_darkmode" align=middle width=29.1030036pt height=22.4657235pt/> points toward node <img src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.6632257pt height=21.6830097pt/> perpendicular to <img src="./tex/2dfa893f7bf3c5c01ceb4788da782da0.svg?invert_in_darkmode" align=middle width=13.3152525pt height=14.6118786pt/> and
likewise <img src="./tex/e7299099834dc07da267043c72b8955e.svg?invert_in_darkmode" align=middle width=30.556614pt height=22.4657235pt/> points toward node <img src="./tex/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710417pt height=21.6830097pt/> perpendicular to <img src="./tex/6ac64257e8429516f022faf3750898f9.svg?invert_in_darkmode" align=middle width=14.7688629pt height=14.6118786pt/>. Call the angle
formed between these two vectors <img src="./tex/fd7c3eb33987d32e5867ff7bb1cacb55.svg?invert_in_darkmode" align=middle width=8.17352745pt height=22.8310566pt/>. So we may write:

<p align="center"><img src="./tex/dde8778be23b02252fad310beca5f28a.svg?invert_in_darkmode" align=middle width=357.5762223pt height=34.7253258pt/></p>


Now notice that the angle between <img src="./tex/2dfa893f7bf3c5c01ceb4788da782da0.svg?invert_in_darkmode" align=middle width=13.3152525pt height=14.6118786pt/> and <img src="./tex/6ac64257e8429516f022faf3750898f9.svg?invert_in_darkmode" align=middle width=14.7688629pt height=14.6118786pt/>, call it <img src="./tex/5e025aa9bc7c6818199db0d8fb82bd3e.svg?invert_in_darkmode" align=middle width=21.27105585pt height=14.1552444pt/>,
is <img src="./tex/01f40bdda54d9ab256ddd5e9a04e371b.svg?invert_in_darkmode" align=middle width=38.22480585pt height=22.8310566pt/>, but more importantly that:
<p align="center"><img src="./tex/5509cfc59c55cd5354fd80489e2d1a13.svg?invert_in_darkmode" align=middle width=234.974784pt height=17.0319402pt/></p>

So, we can rewrite equation the cosine law equation above into:
<p align="center"><img src="./tex/151577143ba9b968bc58f486847c5246.svg?invert_in_darkmode" align=middle width=137.4281601pt height=34.7253258pt/></p>

Now, apply the definition of sine for right triangles:
<p align="center"><img src="./tex/ad301946877f66beddd043039f1845d9.svg?invert_in_darkmode" align=middle width=163.9118646pt height=38.5152603pt/></p>

where <img src="./tex/ddd3bc35b936d6a00e6a81cab0061f32.svg?invert_in_darkmode" align=middle width=14.1220134pt height=22.8310566pt/> is the height of the triangle treating <img src="./tex/2dfa893f7bf3c5c01ceb4788da782da0.svg?invert_in_darkmode" align=middle width=13.3152525pt height=14.6118786pt/> as base, and
likewise for <img src="./tex/6d22be1359e204374e6f0b45e318d561.svg?invert_in_darkmode" align=middle width=15.5756238pt height=22.8310566pt/>. Rewriting the equation above, replacing one of the edge norms,
e.g.\ <img src="./tex/d0cc021c5b0e6847b9929fd5c6ecd826.svg?invert_in_darkmode" align=middle width=30.5755692pt height=24.657534pt/>:
<p align="center"><img src="./tex/99e92f13a02d606f5ce1fd991794f828.svg?invert_in_darkmode" align=middle width=148.3049766pt height=42.9260634pt/></p>


Combine the cosine and sine terms:
<p align="center"><img src="./tex/bd842203a257378897dd2d2499081e97.svg?invert_in_darkmode" align=middle width=127.3093008pt height=34.7253258pt/></p>


Finally, since <img src="./tex/940aebca1f028d93419529cff94a740a.svg?invert_in_darkmode" align=middle width=90.8923323pt height=24.657534pt/>, our constant dot product of these
gradients in our triangle is:
<p align="center"><img src="./tex/40cb87da7035dbcd178fcf1363089dfc.svg?invert_in_darkmode" align=middle width=163.14031305pt height=32.50746015pt/></p>


Similarly, inside the other triangle <img src="./tex/e0563ab36735bded51644cf89825f1d2.svg?invert_in_darkmode" align=middle width=17.6513436pt height=22.4657235pt/> incident
on nodes <img src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.6632257pt height=21.6830097pt/> and <img src="./tex/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710417pt height=21.6830097pt/> with angle <img src="./tex/45f257d25bb292b0f5b8a984f4efbe84.svg?invert_in_darkmode" align=middle width=20.05337235pt height=22.8310566pt/> we have a constant dot
product:
<p align="center"><img src="./tex/4f16a031aba6bcb3e69a09759c051f29.svg?invert_in_darkmode" align=middle width=161.9226114pt height=33.8120871pt/></p>

where <img src="./tex/61e84f854bc6258d4108d08d4c4a0852.svg?invert_in_darkmode" align=middle width=13.2934098pt height=22.4657235pt/> is the area <img src="./tex/e0563ab36735bded51644cf89825f1d2.svg?invert_in_darkmode" align=middle width=17.6513436pt height=22.4657235pt/>.

Recall that <img src="./tex/965b9779f7d35773b987d1080ab428e9.svg?invert_in_darkmode" align=middle width=15.40433235pt height=14.1552444pt/> and <img src="./tex/5e5d7941135a825fa8ab551e0eb76a4f.svg?invert_in_darkmode" align=middle width=16.85794275pt height=14.1552444pt/> are only both nonzero inside these two
triangles, <img src="./tex/04fdc9ebe2a4f15ffdd732d0abd189e1.svg?invert_in_darkmode" align=middle width=18.1520988pt height=22.4657235pt/> and <img src="./tex/e0563ab36735bded51644cf89825f1d2.svg?invert_in_darkmode" align=middle width=17.6513436pt height=22.4657235pt/>.  So, since these constants are inside an
integral over area  we may write:
<p align="center"><img src="./tex/9f8210836fb3ca3da88e84938171d18f.svg?invert_in_darkmode" align=middle width=560.18018265pt height=47.1647583pt/></p>


## Mass matrix

Treated as an _operator_ (i.e., when used multiplied against a vector <img src="./tex/5435f08c809898be0b38d31c7232e2af.svg?invert_in_darkmode" align=middle width=21.8720271pt height=22.5570873pt/>),
the Laplacian matrix <img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979855pt height=22.5570873pt/> computes the local integral of the Laplacian of a
function <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/>. In the energy-based formulation of the smoothing problem this is
not an issue. If we used a similar FEM derivation for the _data term_ we would
get another sparse matrix <img src="./tex/e81a1a11e934fd54074805470b4f9ad1.svg?invert_in_darkmode" align=middle width=76.4345043pt height=26.1773094pt/>:

<p align="center"><img src="./tex/e9b8736a9aa66b1921af8ad5a22571df.svg?invert_in_darkmode" align=middle width=591.3941253pt height=37.3519608pt/></p>

where <img src="./tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg?invert_in_darkmode" align=middle width=17.9451195pt height=22.5570873pt/> as an operator computes the local integral of a function's value
(i.e., <img src="./tex/5a6a7ddb47dbe4427df2767dd2254d3c.svg?invert_in_darkmode" align=middle width=28.44734805pt height=22.5570873pt/>).

This matrix <img src="./tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg?invert_in_darkmode" align=middle width=17.9451195pt height=22.5570873pt/> is often _diagonalized_ or _lumped_ into a diagonal matrix,
even in the context of FEM. So often we will simply set:

<p align="center"><img src="./tex/48f326d21a7a528a0f1aeb6a90755416.svg?invert_in_darkmode" align=middle width=497.2228635pt height=71.0142642pt/></p>

for a mesh with <img src="./tex/0e51a2dede42189d77627c4d742822c3.svg?invert_in_darkmode" align=middle width=14.4331011pt height=14.1552444pt/> triangles.

If we start directly with the continuous smoothing iteration equation, then we
have a point-wise equality. To fit in our integrated Laplacian <img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979855pt height=22.5570873pt/> we should
convert it to a point-wise quantity. From a units perspective, we need to
divide by the local area. This would result in a discrete smoothing iteration
equation:

<p align="center"><img src="./tex/5d53591a59f665c794e0a2f859142b53.svg?invert_in_darkmode" align=middle width=172.30560435pt height=18.3123831pt/></p>


where <img src="./tex/d4fecafea88cb244ca98be86cb2bb526.svg?invert_in_darkmode" align=middle width=65.6583081pt height=26.1773094pt/> is the identity matrix. This equation is _correct_ but
the resulting matrix <img src="./tex/4bdc180d4843f731cfac9a4b56ba3366.svg?invert_in_darkmode" align=middle width=124.58859435pt height=26.7617526pt/> is not symmetric and thus slower
to solve against.

Instead, we could take the healthier view of requiring our smoothing iteration
equation to hold in a locally integrated sense. In this case, we replace mass
matrices on either side:

<p align="center"><img src="./tex/58bea700f9fd133dfde44969079927c0.svg?invert_in_darkmode" align=middle width=165.4333428pt height=18.3123831pt/></p>


Now the system matrix <img src="./tex/c317e35fe20ee616ed6678d04169794b.svg?invert_in_darkmode" align=middle width=99.7712133pt height=22.8310566pt/> will be symmetric and we can use
[Cholesky factorization](https://en.wikipedia.org/wiki/Cholesky_decomposition)
to solve with it.

### Laplace Operator is Intrinsic

The discrete Laplacian operator and its accompanying mass matrix are
_intrinsic_ operators in the sense that they _only_ depend on lengths. In
practical terms, this means we do not need to know _where_ vertices are
actually positioned in space (i.e., <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.55473745pt height=22.5570873pt/>). Rather we only need to know the
relative distances between neighboring vertices (i.e., edge lengths). We do not
even need to know which dimension this mesh is [living
in](https://en.wikipedia.org/wiki/Embedding).

This also means that applying a transformation to a shape that does not change
any lengths on the surface (e.g., bending a sheet of paper) will have no affect
on the Laplacian.

### Data denoising

For the data denoising application, our geometry of the domain is not changing
only the scalar function living upon it. We can build our discrete Laplacian
<img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979855pt height=22.5570873pt/> and mass matrix <img src="./tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg?invert_in_darkmode" align=middle width=17.9451195pt height=22.5570873pt/> and apply the above formula with a chosen <img src="./tex/1b109d8b4484cf614f27126d788c510e.svg?invert_in_darkmode" align=middle width=9.58908225pt height=22.8310566pt/>
parameter.

### Geometric smoothing

For geometric smoothing, the Laplacian operator (both <img src="./tex/545994479b8eb0934d80a97627e63f82.svg?invert_in_darkmode" align=middle width=13.69867125pt height=22.4657235pt/> in the continuous
setting and <img src="./tex/a72b56a59ffa2be52208086ccb6f7562.svg?invert_in_darkmode" align=middle width=36.6208029pt height=22.5570873pt/> in the discrete setting) depend on the geometry of the
surface <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.5022269pt height=22.5570873pt/>. So if the signal <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/> is replaced with the positions of points on
the surface (say, <img src="./tex/8fc2ed1e8ffd9008ef2661b89ca5b616.svg?invert_in_darkmode" align=middle width=71.47064265pt height=26.7617526pt/> in the discrete case), then the smoothing
iteration update rule is a _non-linear_ function if we write it as:

<p align="center"><img src="./tex/c6f234133df83b61064c6cdda4c75eb8.svg?invert_in_darkmode" align=middle width=240.8331519pt height=18.3123831pt/></p>


However, if we assume that small changes in <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.55473745pt height=22.5570873pt/> have a negligible effect on
<img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979855pt height=22.5570873pt/> and <img src="./tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg?invert_in_darkmode" align=middle width=17.9451195pt height=22.5570873pt/> then we can discretize _explicitly_ by computing <img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979855pt height=22.5570873pt/> and <img src="./tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg?invert_in_darkmode" align=middle width=17.9451195pt height=22.5570873pt/>
_before_ performing the update:

<p align="center"><img src="./tex/8c168fd96735d211e2c39766d370606a.svg?invert_in_darkmode" align=middle width=190.9014195pt height=18.3123831pt/></p>


### Why did my mesh disappear?

Repeated application of geometric smoothing may cause the mesh to "disappear".
Actually the updated vertex values are being set to
[NaNs](https://en.wikipedia.org/wiki/NaN) due to degenerate numerics. We are
rebuilding the discrete Laplacian at every new iteration, regardless of the
"quality" of the mesh's triangles. In particular, if a triangle tends to become
skinnier and skinnier during smoothing, what will happen to the cotangents of
its angles?

In "Can Mean-Curvature Flow Be Made Non-Singular?", Kazhdan et al. derive a new
type of geometric flow that is stable (so long as the mesh at time <img src="./tex/1c899e1c767eb4eac89facb5d1f2cb0d.svg?invert_in_darkmode" align=middle width=36.0729369pt height=21.1872144pt/> is
reasonable). Their change is remarkably simple: do not update <img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979855pt height=22.5570873pt/>, only update
<img src="./tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg?invert_in_darkmode" align=middle width=17.9451195pt height=22.5570873pt/>.

## Tasks

### Learn an alternative derivation of cotangent Laplacian

The "cotangent Laplacian" by far the most important tool in geometry
processing. It comes up everywhere. It is important to understand where it
comes from and be able to derive it (in one way or another).

The background section above contains a FEM derivation of the discrete
"cotangnet Laplacian". For this (unmarked) task, read and understand one of the
_other_ derivations listed above.

> **Hint:** The finite-volume method used in [Botsch et al. 2010] is perhaps
> the most accessible alternative.

### White list

  - `igl::doublearea`
  - `igl::edge_lengths`

### Black list

  - `igl::cotmatrix_entries`
  - `igl::cotmatrix`
  - `igl::massmatrix`
  - Trig functions `sin`, `cos`, `tan` etc. (e.g., from `#include <cmath>`)
    _See background notes about "intrinisic"-ness_

### `src/cotmatrix.cpp`

Construct the "cotangent Laplacian" for a mesh with edge lengths `l`. Each
entry in the output sparse, symmetric matrix `L` is given by:

<p align="center"><img src="./tex/b135eeed3e8770ac50d05128df513462.svg?invert_in_darkmode" align=middle width=325.8550152pt height=69.0417981pt/></p>


> Hint: Review the [law of sines](https://en.wikipedia.org/wiki/Law_of_sines)
> and [law of cosines](https://en.wikipedia.org/wiki/Law_of_cosines) and
> [Heron's ancient formula](https://en.wikipedia.org/wiki/Heron's_formula) to
> derive a formula for the cotangent of each triangle angle that _only_ uses
> edge lengths.

### `src/massmatrix.cpp`

Construct the diagonal(ized) mass matrix `M` for a mesh with given face indices
in `F` and edge lengths `l`.

### `src/smooth.cpp`

Given a mesh (`V`,`F`) and data specified per-vertex (`G`), smooth this data
using a single implicit Laplacian smoothing step.

This data could be a scalar field on the surface and smoothing corresponds to
data denoising.

![](images/beetle-data-denoising.gif)

Or the data could be the vector field of the surface's own geometry. This
corresponds to geometric smoothing.

![](images/sphere-geometric-smoothing.gif)

