# Geometry Processing – Smoothing

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
differential equation we will start with sets the change in signal value $u$
over time [proportional
to](https://en.wikipedia.org/wiki/Proportionality_(mathematics)) the
[Laplacian](https://en.wikipedia.org/wiki/Laplace_operator) of the signal $∆u$
(for now, roughly the second derivative of the signal as we move _on_ the
surface):

\\[
\frac{∂ u}{∂ t} = λ ∆ u,
\\]

where the scalar parameter $λ$ controls the rate of smoothing.

When the signal is the surface geometry, we call this a [geometric
flow](https://en.wikipedia.org/wiki/Geometric_flow).

There are various ways to motivate this choice of flow for
data-/geometry-smoothing. Let us consider one way that will introduce the
Laplacian as a form of local averaging.

Given a noisy signal $f$, intuitively we can _smooth_ $f$ by averaging every
value with its neighbors' values. In continuous math, we might write that the
smoothed value $u(\x)$ at any point on our surface $\x ∈ \S$ should be equal to
the average value of some small
[ball](https://en.wikipedia.org/wiki/Ball_(mathematics)) of nearby points:

\\[
u(\x) = \frac{1}{|B(\x))|} ∫_{B(\x)} f(\z) \;d\z,
\\]

If the ball $B(\x)$ is small, then we will have to repeat this averaging many
times to see a global smoothing effect. Hence, we can write that the current
value $u^t$ _flows_ toward smooth solution by small steps $δt$ in time:

\\[
u^{t+δt}(\x) = \frac{1}{|B(\x))|} ∫_{B(\x)} u^t(\z) \;d\z.
\\]

Subtracting the current value $u^t(\x)$ from both sides and introducing a
flow-speed parameter $λ$ we have a flow equation
describing the change in value as an integral of relative values:

\\[
\frac{∂ u}{∂ t} 
  = λ \frac{1}{|B(\x))|} ∫_{B(\x)} u(\z)-u(\x) \;d\z.
\\]

For harmonic functions, $∆u =0$, this integral becomes zero in the limit as the
radius of the ball shrinks to zero via satisfaction of the
[mean value
theorem](https://en.wikipedia.org/wiki/Harmonic_function#The_mean_value_property).
It follows for a non-harmonic $∆u ≠ 0$ this integral is equal to the Laplacian
of the $u$, so we have arrived at our flow equation:

\\[
\frac{∂ u}{∂ t} 
  = \lim_{|B(\x)| → 0}  λ \frac{1}{|B(\x))|} ∫_{B(\x)} u(\z)-u(\x) \;d\z 
  = λ ∆ u.
\\]

### Energy-based formulation

Alternatively, we can think of a single smoothing operation as the solution to
an energy minimization problem. If $f$ is our noisy signal over the surface,
then we want to find a signal $u$ such that it simultaneously minimizes its
difference with $f$ and minimizes its variation over the surface:

\\[
u^* 
  = \argmin_u E_(u) 
  = 
  \argmin_u ½∫_\S \underbrace{(f-u)²}_\text{data term} + 
  λ\underbrace{‖∇u‖²}_\text{smoothness term} \;dA,
\\]

where again the scalar parameter $λ$ controls the rate of smoothing. This
energy-based formulation is equivalent to the flow-based formulation.
Minimizing these energies is identical to stepping forward one temporal unit in
the flow.

#### Calculus of variations
In the smooth setting, our energy $E$ is a function that measures scalar value
of a given function _u_, making it a
[functional](https://en.wikipedia.org/wiki/Functional_(mathematics)). To
understand how to _minimize_ a functional with respect to an unknown function,
we will need concepts from the [calculus of
variations](https://en.wikipedia.org/wiki/Calculus_of_variations).

We are used to working with minimizing quadratic _functions_ with respect to a
discrete set of variables, where the minimum is obtained when the gradient of
the energy with respect to the variables is zero.

In our case, the functional $E(u)$ is quadratic in $u$ (recall that the
[gradient operator](https://en.wikipedia.org/wiki/Gradient) $∇$ is a linear
operator). The function $u$ that minimizes $E(u)$ will be obtained when any
small change or _variation_ in $u$ has no change on the energy values. To
create a small change in a function $u$ we will add another function $v$ times
a infinitesimal scalar $ε$. If $E(u)$ is minimized for a function $w$ and we
are given another arbitrary function $v$, then let us define a function new
function 

\\[
Φ(ε) = E(w+εv) = ½∫_\S (f-w+εv)² + λ ‖∇w + ε∇v‖²  \;dA,
\\]
where we observe that $Φ$ is quadratic in $ε$.

Since $E(w)$ is minimal then $Φ$ is minimized when $ε$ is zero, and if $Φ$ is
minimal at $ε=0$, then the derivative of $Φ$ with respect $ε$ must be zero:

\\[
\begin{align}
0 & = \left.\frac{∂Φ}{∂ε} \right|_{ε = 0},\\
  & = \left.\frac{∂}{∂ε} ½∫_\S (f-w-εv)² + λ ‖∇w + ε∇v‖²\;dA, \right|_{ε = 0} \\
  & = \left.\frac{∂}{∂ε} ½∫_\S
    f^2 - 2wf - 2εfv +w²+2εvw +ε²v² + λ ‖∇w‖² + λ2ε∇v⋅∇w + λ ε²‖∇w‖² \;dA \right|_{ε = 0}\\
  & = \left.∫_\S -fv + vw +2εvw  + λ∇v⋅∇w + λ ε‖∇w‖² \;dA \right|_{ε = 0}\\
  & = ∫_\S v(w-f)  + λ∇v⋅∇w \;dA.
\end{align}
\\]

The choice of "test" function $v$ was arbitrary, so this must hold for any
(reasonable) choice of $v$:

\\[
0 = ∫_\S v(w-f)  + λ∇v⋅∇w \;dA \quad ∀ v.
\\]

It is difficult to claim much about $w$ from this equation directly because
derivatives of $v$ are still involved. We can _move_ a derivative from $v$ to a
$w$ by applying [Green's first
identity](https://en.wikipedia.org/wiki/Green's_identities):

\\[
0 = ∫_\S v(w-f)  - λv∆w \;dA \quad (+  \text{boundary term} )\quad ∀ v,
\\]
where we choose to _ignore_ the boundary terms (for now) or equivalently we
agree to work on _closed_ surfaces $\S$.

Since this equality must hold of _any_ $v$ let us consider functions that are
little ["blips"](https://en.wikipedia.org/wiki/Bump_function) centered at any
arbitrary point $\x ∈ \S$. A function $v$ that is one at $\x$ and quickly
decays to zero everywhere else. To satisfy the equation above at $\x$ with this
blip $v$ we must have that:

\\[
w(\x)-f(\x) = λ∆w(\x).
\\]
The choice of $\x$ was arbitrary so this must hold _everywhere_.

Because we invoke _variations_ to arrive at this equation, we call the
_energy-based_ formulation a _variational formulation_.

### Implicit smoothing iteration

Now we understand that the flow-based formulation and the variational
formulation lead to the same system, let us concretely write out the implicit
smoothing step. 

Letting $u^0 = f$ we compute a new smoothed function $u^{t+1}$ given the
current solution $u^t$ by solving the _linear_ system of equations:

\\[
u^t(\x) = (\text{id}-λ∆)u^{t+1}(\x), \quad ∀ \x ∈ \S
\\]
where $\text{id}$ is the [identity
operator](https://en.wikipedia.org/wiki/Identity_function). In the discrete
case, we will need discrete approximations of the $\text{id}$ and $∆$
operators.

## Discrete Laplacian

There are many ways to derive a discrete
approximation of the Laplacian $∆$ operator on a triangle mesh using:

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
 - [gradient of surface area → mean curvature
   flow](https://en.wikipedia.org/wiki/Mean_curvature_flow)
    - "Computing Discrete Minimal Surfaces and Their Conjugates" [Pinkall &
      Polthier 1993]

All of these techniques will produce the _same_ sparse _Laplacian matrix_ $\L ∈ \R^{n×n}$
for a mesh with $n$ vertices. 

### Finite volume derivation of the discrete Laplacian

With minimal background theory, we can use the [finite volume
method](https://en.wikipedia.org/wiki/Finite_volume_method) to define our
discrete approximation of the Laplacian operator $∆$ on a given triangle mesh.

The rough idea of the finite volume method is to approximate a smooth function
over a triangle mesh by considering integral averages in small regions around
each mesh vertex.

We want to approximate the Laplacian of a function $∆u$. Let us consider $u$ to
be [piecewise-linear](https://en.wikipedia.org/wiki/Piecewise_linear_function)
represented by scalar values at each vertex, collected in $\u ∈ \R^n$.

We will approximate $∆u$ as the [integral
average](https://en.wikipedia.org/wiki/Mean_of_a_function) of $∆u$ in the local
neighborhood $N(\v_i)$ of each mesh vertex:

\\[
∆u(\v_i) ≈ \frac{1}{\text{Area}(N(\v_i))} ∫_{N(\v)} ∆u \;dA.
\\]

**Voronoi neighborhood:** One of the simplest choices for _dividing up_ our
mesh into  regions around each vertex is to consider the [Voronoi
diagram](https://en.wikipedia.org/wiki/Voronoi_diagram). The neighborhood
$N(\v_i)$ around a vertex $\v_i$ is defined as all of the points that are
closer to $\v_i$ than any other vertex $\v_j$.

> For triangle mesh surfaces in $\R³$ the Voronoi diagram has to be defined
> with respect to [distance _on_ the
> mesh](https://en.wikipedia.org/wiki/Geodesic) rather than usual distance in
> $\R³$. But luckily we will not have to explicitly compute these distances.



## Data denoising

Let us start with the data denoising application. In this case, our geometry of
the domain is not changing. 

### Why did my mesh disappear?

## Cotangent Laplacian

  - Finite Volume (Mark Meyer/PMP book)
  - FEM (Courant etc.)
  - Mean curvature normal, Area gradient (Desbrun/Pinkall)
  - DEC Desbrun/Hirani
  - GBC Wachspress/Floater

## Tasks

### `src/cotmatrix.cpp`

> Hint: Review the [law of sines](https://en.wikipedia.org/wiki/Law_of_sines)
> and [law of cosines](https://en.wikipedia.org/wiki/Law_of_cosines) and
> [Heron's ancient formula](https://en.wikipedia.org/wiki/Heron's_formula) to
> derive a formula for the cotangent of each triangle angle that _only_ uses
> edge lengths.

### White list

  - `igl::doublearea`
  - `igl::edge_lengths`

### Black list


  - `igl::cotmatrix_entries`
  - `igl::cotmatrix`
  - `igl::massmatrix`
  - Trig functions `sin`, `cos`, `tan` etc. (e.g., from `#include <cmath>`)
