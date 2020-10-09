# Geometry Processing Smoothing

> **To get started:** Fork this repository then issue
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
differential equation we will start with sets the change in signal value <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/>
over time [proportional
to](https://en.wikipedia.org/wiki/Proportionality_(mathematics)) the
[Laplacian](https://en.wikipedia.org/wiki/Laplace_operator) of the signal <img alt="$\Delta u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/520957eebcbc4341addf55ec1a6ef0ab.svg" align="middle" width="23.10894464999999pt" height="22.465723500000017pt"/>
(for now, roughly the second derivative of the signal as we move _on_ the
surface):

<p align="center"><img alt="$$&#10;\frac{\partial  u}{\partial  t} = {\lambda} \Delta  u,&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/eb0a14f1824b0e6901ad51c2b68e8ee7.svg" align="middle" width="80.20514865pt" height="33.81208709999999pt"/></p>


where the scalar parameter <img alt="${\lambda}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/1b109d8b4484cf614f27126d788c510e.svg" align="middle" width="9.58908224999999pt" height="22.831056599999986pt"/> controls the rate of smoothing.

When the signal is the surface geometry, we call this a [geometric
flow](https://en.wikipedia.org/wiki/Geometric_flow).

There are various ways to motivate this choice of flow for
data-/geometry-smoothing. Let us consider one way that will introduce the
Laplacian as a form of local averaging.

Given a noisy signal <img alt="$f$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/190083ef7a1625fbc75f243cffb9c96d.svg" align="middle" width="9.81741584999999pt" height="22.831056599999986pt"/>, intuitively we can _smooth_ <img alt="$f$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/190083ef7a1625fbc75f243cffb9c96d.svg" align="middle" width="9.81741584999999pt" height="22.831056599999986pt"/> by averaging every
value with its neighbors' values. In continuous math, we might write that the
smoothed value <img alt="$u(\mathbf{x})$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/193e99dd587ecc23cdd8853e4cd46a76.svg" align="middle" width="32.17282199999999pt" height="24.65753399999998pt"/> at any point on our surface <img alt="$\mathbf{x} \in  \mathbf{S}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/191fb656c7f4f20997cf240f7af26caf.svg" align="middle" width="40.57048049999999pt" height="22.55708729999998pt"/> should be equal to
the average value of some small
[ball](https://en.wikipedia.org/wiki/Ball_(mathematics)) of nearby points:

<p align="center"><img alt="$$&#10;u(\mathbf{x}) = \frac{1}{|B(\mathbf{x}))|} \int _{B(\mathbf{x})} f(\mathbf{z}) \;d\mathbf{z},&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/a69ec5fa9866f2be8caf1efd6fdc8d11.svg" align="middle" width="210.76788645pt" height="40.2286731pt"/></p>


If the ball <img alt="$B(\mathbf{x})$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/06e8e88b4ac31efcc03ecf253f9609ba.svg" align="middle" width="36.05595344999999pt" height="24.65753399999998pt"/> is small, then we will have to repeat this averaging many
times to see a global smoothing effect. Hence, we can write that the current
value <img alt="$u^t$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2781fb506c10958a89548aa070f15594.svg" align="middle" width="14.37606554999999pt" height="26.085962100000025pt"/> _flows_ toward smooth solution by small steps <img alt="${\delta}t$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/c6bcca8f259a226ca0748027158124bd.svg" align="middle" width="13.864169549999989pt" height="22.831056599999986pt"/> in time:

<p align="center"><img alt="$$&#10;u^{t+{\delta}t}(\mathbf{x}) = \frac{1}{|B(\mathbf{x}))|} \int _{B(\mathbf{x})} u^t(\mathbf{z}) \;d\mathbf{z}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e8b064f3dd720aa90964378e9d8b8f47.svg" align="middle" width="243.3853983pt" height="40.2286731pt"/></p>


Subtracting the current value <img alt="$u^t(\mathbf{x})$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/4813491b0a749868771cab60ce88649f.svg" align="middle" width="37.960510499999984pt" height="26.085962100000025pt"/> from both sides and introducing a
flow-speed parameter <img alt="${\lambda}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/1b109d8b4484cf614f27126d788c510e.svg" align="middle" width="9.58908224999999pt" height="22.831056599999986pt"/> we have a flow equation
describing the change in value as an integral of relative values:

<p align="center"><img alt="$$&#10;\frac{\partial  u}{\partial  t} &#10;  = {\lambda} \frac{1}{|B(\mathbf{x}))|} \int _{B(\mathbf{x})} (u(\mathbf{z})-u(\mathbf{x})) \;d\mathbf{z}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/32bab51d0ad9d94ce1fe47d67e452dbb.svg" align="middle" width="271.11006119999996pt" height="40.3902642pt"/></p>


For harmonic functions, <img alt="$\Delta u =0$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/bb23bb09b67bb167730c76ff706043de.svg" align="middle" width="53.245785449999985pt" height="22.465723500000017pt"/>, this integral becomes zero in the limit as the
radius of the ball shrinks to zero via satisfaction of the
[mean value
theorem](https://en.wikipedia.org/wiki/Harmonic_function#The_mean_value_property).
It follows for a non-harmonic <img alt="$\Delta u \ne  0$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/74a93b9b98af7c7b321f354c4fcf8ea2.svg" align="middle" width="53.245785449999985pt" height="22.831056599999986pt"/> this integral is equal to the Laplacian
of the <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/>, so we have arrived at our flow equation:

<p align="center"><img alt="$$&#10;\frac{\partial  u}{\partial  t} &#10;  = \lim_{|B(\mathbf{x})| \Rightarrow  0}  {\lambda} \frac{1}{|B(\mathbf{x}))|} \int _{B(\mathbf{x})} (u(\mathbf{z})-u(\mathbf{x})) \;d\mathbf{z} &#10;  = {\lambda} \Delta  u.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/fd1ea805db98db0bfe3a05020cbc7771.svg" align="middle" width="384.55396979999995pt" height="40.3902642pt"/></p>


### Energy-based formulation

Alternatively, we can think of a single smoothing operation as the solution to
an energy minimization problem. If <img alt="$f$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/190083ef7a1625fbc75f243cffb9c96d.svg" align="middle" width="9.81741584999999pt" height="22.831056599999986pt"/> is our noisy signal over the surface,
then we want to find a signal <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/> such that it simultaneously minimizes its
difference with <img alt="$f$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/190083ef7a1625fbc75f243cffb9c96d.svg" align="middle" width="9.81741584999999pt" height="22.831056599999986pt"/> and minimizes its variation over the surface:

<p align="center"><img alt="$$&#10;u^* &#10;  = \mathop{\text{argmin}}_u E_(u) &#10;  = &#10;  \mathop{\text{argmin}}_u \frac12 \int _\mathbf{S} ( \underbrace{(f-u)^{2}}_\text{data} + &#10;  \underbrace{{\lambda}\| {\nabla}u\| ^{2}}_\text{smoothness} )\;dA,&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2103cd36ae87c74dbd5cf1297f2bb1da.svg" align="middle" width="419.77996004999994pt" height="49.05031724999999pt"/></p>


where again the scalar parameter <img alt="${\lambda}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/1b109d8b4484cf614f27126d788c510e.svg" align="middle" width="9.58908224999999pt" height="22.831056599999986pt"/> controls the rate of smoothing. This
energy-based formulation is equivalent to the flow-based formulation.
Minimizing these energies is identical to stepping forward one temporal unit in
the flow.

#### Calculus of variations
In the smooth setting, our energy <img alt="$E$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/84df98c65d88c6adf15d4645ffa25e47.svg" align="middle" width="13.08219659999999pt" height="22.465723500000017pt"/> is a function that measures scalar value
of a given function _u_, making it a
[functional](https://en.wikipedia.org/wiki/Functional_(mathematics)). To
understand how to _minimize_ a functional with respect to an unknown function,
we will need concepts from the [calculus of
variations](https://en.wikipedia.org/wiki/Calculus_of_variations).

We are used to working with minimizing quadratic _functions_ with respect to a
discrete set of variables, where the minimum is obtained when the gradient of
the energy with respect to the variables is zero.

In our case, the functional <img alt="$E(u)$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/36f9bbd82b66589fd9772b7889d95932.svg" align="middle" width="35.27788274999999pt" height="24.65753399999998pt"/> is quadratic in <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/> (recall that the
[gradient operator](https://en.wikipedia.org/wiki/Gradient) <img alt="${\nabla}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/289bc3c68d369342d98da6ff43ec5b9f.svg" align="middle" width="13.69867124999999pt" height="22.465723500000017pt"/> is a linear
operator). The function <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/> that minimizes <img alt="$E(u)$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/36f9bbd82b66589fd9772b7889d95932.svg" align="middle" width="35.27788274999999pt" height="24.65753399999998pt"/> will be obtained when any
small change or _variation_ in <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/> has no change on the energy values. To
create a small change in a function <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/> we will add another function <img alt="$v$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/> times
a infinitesimal scalar <img alt="${\epsilon}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/7a676caa07476cecc1fef337bd6a860e.svg" align="middle" width="6.672392099999992pt" height="14.15524440000002pt"/>. If <img alt="$E(u)$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/36f9bbd82b66589fd9772b7889d95932.svg" align="middle" width="35.27788274999999pt" height="24.65753399999998pt"/> is minimized for a function <img alt="$w$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/31fae8b8b78ebe01cbfbe2fe53832624.svg" align="middle" width="12.210846449999991pt" height="14.15524440000002pt"/> and we
are given another arbitrary function <img alt="$v$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/>, then let us define a function new
function 

<p align="center"><img alt="$$&#10;{\phi}({\epsilon}) = E(w+{\epsilon}v) = \frac12 \int _\mathbf{S} ((f-w+{\epsilon}v)^{2} + {\lambda} \| {\nabla}w + {\epsilon}{\nabla}v\| ^{2})  \;dA,&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/1676aad3b697e1ce24e6fdf9e3b80555.svg" align="middle" width="448.4368152pt" height="37.3519608pt"/></p>

where we observe that <img alt="${\phi}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/eb62325152b7dd734820edcdd33d2950.svg" align="middle" width="9.794543549999991pt" height="22.831056599999986pt"/> is quadratic in <img alt="${\epsilon}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/7a676caa07476cecc1fef337bd6a860e.svg" align="middle" width="6.672392099999992pt" height="14.15524440000002pt"/>.

Since <img alt="$E(w)$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/21d2cb0fed44339d6eb9fd1c6788522f.svg" align="middle" width="38.078455799999986pt" height="24.65753399999998pt"/> is minimal then <img alt="${\phi}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/eb62325152b7dd734820edcdd33d2950.svg" align="middle" width="9.794543549999991pt" height="22.831056599999986pt"/> is minimized when <img alt="${\epsilon}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/7a676caa07476cecc1fef337bd6a860e.svg" align="middle" width="6.672392099999992pt" height="14.15524440000002pt"/> is zero, and if <img alt="${\phi}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/eb62325152b7dd734820edcdd33d2950.svg" align="middle" width="9.794543549999991pt" height="22.831056599999986pt"/> is
minimal at <img alt="${\epsilon}=0$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/00e14c64b601ac983b4364aa60fa1dc6.svg" align="middle" width="36.80923124999999pt" height="21.18721440000001pt"/>, then the derivative of <img alt="${\phi}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/eb62325152b7dd734820edcdd33d2950.svg" align="middle" width="9.794543549999991pt" height="22.831056599999986pt"/> with respect <img alt="${\epsilon}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/7a676caa07476cecc1fef337bd6a860e.svg" align="middle" width="6.672392099999992pt" height="14.15524440000002pt"/> must be zero:


The choice of "test" function <img alt="$v$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/> was arbitrary, so this must hold for any
(reasonable) choice of <img alt="$v$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/>:

<p align="center"><img alt="$$&#10;0 = \int _\mathbf{S} (v(w-f)  + {\lambda}{\nabla}v\cdot {\nabla}w) \;dA \quad \forall  v.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/f918be183bda0a2a6d22fe2690f24395.svg" align="middle" width="278.4682362pt" height="37.3519608pt"/></p>


It is difficult to claim much about <img alt="$w$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/31fae8b8b78ebe01cbfbe2fe53832624.svg" align="middle" width="12.210846449999991pt" height="14.15524440000002pt"/> from this equation directly because
derivatives of <img alt="$v$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/> are still involved. We can _move_ a derivative from <img alt="$v$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/> to a
<img alt="$w$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/31fae8b8b78ebe01cbfbe2fe53832624.svg" align="middle" width="12.210846449999991pt" height="14.15524440000002pt"/> by applying [Green's first
identity](https://en.wikipedia.org/wiki/Green's_identities):

<p align="center"><img alt="$$&#10;0 = \int _\mathbf{S} (v(w-f)  - {\lambda}v\Delta w )\;dA \quad (+  \text{boundary term} )\quad \forall  v,&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/04c2c8b27440c9aea318748e3a1a9ceb.svg" align="middle" width="402.7609938pt" height="37.3519608pt"/></p>

where we choose to _ignore_ the boundary terms (for now) or equivalently we
agree to work on _closed_ surfaces <img alt="$\mathbf{S}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/>.

Since this equality must hold of _any_ <img alt="$v$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/> let us consider functions that are
little ["blips"](https://en.wikipedia.org/wiki/Bump_function) centered at any
arbitrary point <img alt="$\mathbf{x} \in  \mathbf{S}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/191fb656c7f4f20997cf240f7af26caf.svg" align="middle" width="40.57048049999999pt" height="22.55708729999998pt"/>. A function <img alt="$v$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/> that is one at <img alt="$\mathbf{x}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/b0ea07dc5c00127344a1cad40467b8de.svg" align="middle" width="9.97711604999999pt" height="14.611878600000017pt"/> and quickly
decays to zero everywhere else. To satisfy the equation above at <img alt="$\mathbf{x}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/b0ea07dc5c00127344a1cad40467b8de.svg" align="middle" width="9.97711604999999pt" height="14.611878600000017pt"/> with this
blip <img alt="$v$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/> we must have that:

<p align="center"><img alt="$$&#10;w(\mathbf{x})-f(\mathbf{x}) = {\lambda}\Delta w(\mathbf{x}).&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e94a6a9790d776ead8b76d74b10cc90f.svg" align="middle" width="172.38955305pt" height="16.438356pt"/></p>

The choice of <img alt="$\mathbf{x}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/b0ea07dc5c00127344a1cad40467b8de.svg" align="middle" width="9.97711604999999pt" height="14.611878600000017pt"/> was arbitrary so this must hold _everywhere_.

Because we invoke _variations_ to arrive at this equation, we call the
_energy-based_ formulation a _variational formulation_.

### Implicit smoothing iteration

Now we understand that the flow-based formulation and the variational
formulation lead to the same system, let us concretely write out the implicit
smoothing step. 

Letting <img alt="$u^0 = f$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/eb367462f49f6165b7702489cd3ad4ea.svg" align="middle" width="48.51977954999998pt" height="26.76175259999998pt"/> we compute a new smoothed function <img alt="$u^{t+1}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/74341d0452aa58e0f8653b6058b1f3f6.svg" align="middle" width="31.01998514999999pt" height="26.76175259999998pt"/> given the
current solution <img alt="$u^t$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2781fb506c10958a89548aa070f15594.svg" align="middle" width="14.37606554999999pt" height="26.085962100000025pt"/> by solving the _linear_ system of equations:

<p align="center"><img alt="$$&#10;u^t(\mathbf{x}) = (\text{id}-{\lambda}\Delta )u^{t+1}(\mathbf{x}), \quad \forall  \mathbf{x} \in  \mathbf{S}&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/ca2a12218ffab9d11bc2166ed9ad9083.svg" align="middle" width="257.79280725pt" height="18.312383099999998pt"/></p>

where <img alt="$\text{id}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2b41e5149bbd9bd2a0b7170924be9b80.svg" align="middle" width="13.69867124999999pt" height="22.831056599999986pt"/> is the [identity
operator](https://en.wikipedia.org/wiki/Identity_function). In the discrete
case, we will need discrete approximations of the <img alt="$\text{id}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2b41e5149bbd9bd2a0b7170924be9b80.svg" align="middle" width="13.69867124999999pt" height="22.831056599999986pt"/> and <img alt="$\Delta $" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/545994479b8eb0934d80a97627e63f82.svg" align="middle" width="13.69867124999999pt" height="22.465723500000017pt"/>
operators.

## Discrete Laplacian

There are many ways to derive a discrete
approximation of the Laplacian <img alt="$\Delta $" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/545994479b8eb0934d80a97627e63f82.svg" align="middle" width="13.69867124999999pt" height="22.465723500000017pt"/> operator on a triangle mesh using:

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

All of these techniques will produce the _same_ sparse _Laplacian matrix_ <img alt="$\mathbf{L} \in  \mathbb{R}^{n\times n}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/741ed3ef6a246c4afd8317583b95f093.svg" align="middle" width="69.85918335pt" height="26.17730939999998pt"/>
for a mesh with <img alt="$n$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/55a049b8f161ae7cfeb0197d75aff967.svg" align="middle" width="9.86687624999999pt" height="14.15524440000002pt"/> vertices. 

### Finite element derivation of the discrete Laplacian

We want to approximate the Laplacian of a function <img alt="$\Delta u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/520957eebcbc4341addf55ec1a6ef0ab.svg" align="middle" width="23.10894464999999pt" height="22.465723500000017pt"/>. Let us consider <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/> to
be [piecewise-linear](https://en.wikipedia.org/wiki/Piecewise_linear_function)
represented by scalar values at each vertex, collected in <img alt="$\mathbf{u} \in  \mathbb{R}^n$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/1a693d20e239e8c342dd2a059e9770c6.svg" align="middle" width="50.591567399999995pt" height="22.648391699999998pt"/>.

Any piecewise-linear function can be expressed as a sum of values at mesh
vertices times corresponding piecewise-linear basis functions  (a.k.a hat
functions, <img alt="${\varphi}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/965b9779f7d35773b987d1080ab428e9.svg" align="middle" width="15.404332349999988pt" height="14.15524440000002pt"/>):


![](images/hat-function.png)

By plugging this definition into our smoothness energy above, we have discrete
energy that is quadratic in the values at each mesh vertex:



By defining <img alt="${\varphi}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/965b9779f7d35773b987d1080ab428e9.svg" align="middle" width="15.404332349999988pt" height="14.15524440000002pt"/> as piecewise-linear hat functions, the values in the system
matrix <img alt="$L_{ij}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/58e56656e6e782a9941daf494e8471eb.svg" align="middle" width="21.942650399999987pt" height="22.465723500000017pt"/> are uniquely determined by the geometry of the underlying mesh.
These values are famously known as _cotangent weights_. "Cotangent"
because, as we will shortly see, of their trigonometric formulae and "weights"
because as a matrix <img alt="$\mathbf{L}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/80637df1ca7533740cc7b3fdd1ab540b.svg" align="middle" width="11.36979854999999pt" height="22.55708729999998pt"/> they define a weighted [graph
Laplacian](https://en.wikipedia.org/wiki/Laplacian_matrix) for the given mesh.
Graph Laplacians are employed often in geometry processing, and often in
discrete contexts ostensibly disconnected from FEM. The choice or manipulation
of Laplacian weights and subsequent use as a discrete Laplace operator has been
a point of controversy in geometry processing research (see "Discrete laplace
operators: no free lunch" [Wardetzky et al. 2007]).


We first notice that <img alt="${\nabla}{\varphi}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/7d9f97a055843950ec03bc359085c06e.svg" align="middle" width="29.10300359999999pt" height="22.465723500000017pt"/> are constant on each triangle, and only nonzero on
triangles incident on node <img alt="$i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663225699999989pt" height="21.68300969999999pt"/>. For such a triangle, <img alt="$T_{\alpha}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/04fdc9ebe2a4f15ffdd732d0abd189e1.svg" align="middle" width="18.15209879999999pt" height="22.465723500000017pt"/>, this <img alt="${\nabla}{\varphi}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/7d9f97a055843950ec03bc359085c06e.svg" align="middle" width="29.10300359999999pt" height="22.465723500000017pt"/> points
perpendicularly from the opposite edge <img alt="$e_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/b95c2b0aab2482e5bebd25332a4bbde0.svg" align="middle" width="12.30503669999999pt" height="14.15524440000002pt"/> with inverse magnitude equal to
the height <img alt="$h$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2ad9d098b937e46f9f58968551adac57.svg" align="middle" width="9.47111549999999pt" height="22.831056599999986pt"/> of the triangle treating that opposite edge as base:
<p align="center"><img alt="\begin{equation} \|{\nabla}{\varphi}_i\| = \frac{1}{h} = \frac{\|\mathbf{e}_i\|}{2A}, \end{equation}" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/ea89c557e67ed7e586f007c163220bc5.svg" align="middle" width="421.48790805pt" height="34.7253258pt"/></p>
where <img alt="$\mathbf{e}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2dfa893f7bf3c5c01ceb4788da782da0.svg" align="middle" width="13.31525249999999pt" height="14.611878600000017pt"/> is the edge <img alt="$e_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/b95c2b0aab2482e5bebd25332a4bbde0.svg" align="middle" width="12.30503669999999pt" height="14.15524440000002pt"/> as a vector and <img alt="$A$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg" align="middle" width="12.32879834999999pt" height="22.465723500000017pt"/> is the area of the triangle.

![Left: the gradient <img alt="${\nabla} {\varphi}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/461fd2586080bcfe7ac91515a525affe.svg" align="middle" width="29.10300359999999pt" height="22.465723500000017pt"/> of a hat function <img alt="${\varphi}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/965b9779f7d35773b987d1080ab428e9.svg" align="middle" width="15.404332349999988pt" height="14.15524440000002pt"/> is piecewise-constant and
points perpendicular to opposite edges. Right: hat function gradients <img alt="${\nabla} {\varphi}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/461fd2586080bcfe7ac91515a525affe.svg" align="middle" width="29.10300359999999pt" height="22.465723500000017pt"/>
and <img alt="${\nabla} {\varphi}_j$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e7299099834dc07da267043c72b8955e.svg" align="middle" width="30.556613999999993pt" height="22.465723500000017pt"/> of neighboring nodes meet at angle <img alt="${\theta} = {\pi} -&#10;{\alpha}_{ij}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/af94b55c890d07e5e8a5cbccff7647c0.svg" align="middle" width="81.41349149999998pt" height="22.831056599999986pt"/>.](images/hat-function-gradient.png)

Now, consider two neighboring nodes <img alt="$i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663225699999989pt" height="21.68300969999999pt"/> and <img alt="$j$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/36b5afebdba34564d884d347484ac0c7.svg" align="middle" width="7.710416999999989pt" height="21.68300969999999pt"/>, connected by some edge
<img alt="$\mathbf{e}_{ij}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/902505b1f10837bf329ffbeb68762c74.svg" align="middle" width="19.41976244999999pt" height="14.611878600000017pt"/>. Then <img alt="${\nabla}{\varphi}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/7d9f97a055843950ec03bc359085c06e.svg" align="middle" width="29.10300359999999pt" height="22.465723500000017pt"/> points toward node <img alt="$i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663225699999989pt" height="21.68300969999999pt"/> perpendicular to <img alt="$\mathbf{e}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2dfa893f7bf3c5c01ceb4788da782da0.svg" align="middle" width="13.31525249999999pt" height="14.611878600000017pt"/> and
likewise <img alt="${\nabla} {\varphi}_j$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e7299099834dc07da267043c72b8955e.svg" align="middle" width="30.556613999999993pt" height="22.465723500000017pt"/> points toward node <img alt="$j$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/36b5afebdba34564d884d347484ac0c7.svg" align="middle" width="7.710416999999989pt" height="21.68300969999999pt"/> perpendicular to <img alt="$\mathbf{e}_j$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6ac64257e8429516f022faf3750898f9.svg" align="middle" width="14.76886289999999pt" height="14.611878600000017pt"/>. Call the angle
formed between these two vectors <img alt="${\theta}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/fd7c3eb33987d32e5867ff7bb1cacb55.svg" align="middle" width="8.17352744999999pt" height="22.831056599999986pt"/>. So we may write:

<p align="center"><img alt="$$&#10;{\nabla} {\varphi}_i \cdot  {\nabla} {\varphi}_j = \|{\nabla} {\varphi}_i\| \|{\nabla} {\varphi}_j\| \cos {\theta} =&#10;\frac{\|\mathbf{e}_j\|}{2A}\frac{\|\mathbf{e}_i\|}{2A} \cos {\theta}. &#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/dde8778be23b02252fad310beca5f28a.svg" align="middle" width="357.5762223pt" height="34.7253258pt"/></p>


Now notice that the angle between <img alt="$\mathbf{e}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2dfa893f7bf3c5c01ceb4788da782da0.svg" align="middle" width="13.31525249999999pt" height="14.611878600000017pt"/> and <img alt="$\mathbf{e}_j$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6ac64257e8429516f022faf3750898f9.svg" align="middle" width="14.76886289999999pt" height="14.611878600000017pt"/>, call it <img alt="${\alpha}_{ij}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/5e025aa9bc7c6818199db0d8fb82bd3e.svg" align="middle" width="21.27105584999999pt" height="14.15524440000002pt"/>,
is <img alt="${\pi} - {\theta}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/01f40bdda54d9ab256ddd5e9a04e371b.svg" align="middle" width="38.22480584999999pt" height="22.831056599999986pt"/>, but more importantly that:
<p align="center"><img alt="$$&#10;\cos {\theta} = - \cos \left({\pi} - {\theta}\right) = -\cos {\alpha}_{ij}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/5509cfc59c55cd5354fd80489e2d1a13.svg" align="middle" width="234.974784pt" height="17.031940199999998pt"/></p>

So, we can rewrite equation the cosine law equation above into:
<p align="center"><img alt="$$&#10;-\frac{\|\mathbf{e}_j\|}{2A}\frac{\|\mathbf{e}_i\|}{2A} \cos &#10;{\alpha}_{ij}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/151577143ba9b968bc58f486847c5246.svg" align="middle" width="137.42816009999999pt" height="34.7253258pt"/></p>

Now, apply the definition of sine for right triangles:
<p align="center"><img alt="$$&#10;\sin {\alpha}_{ij} = \frac{h_j}{\|\mathbf{e}_i\|} = \frac{h_i}{\|\mathbf{e}_j\|},&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/ad301946877f66beddd043039f1845d9.svg" align="middle" width="163.9118646pt" height="38.5152603pt"/></p>

where <img alt="$h_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/ddd3bc35b936d6a00e6a81cab0061f32.svg" align="middle" width="14.12201339999999pt" height="22.831056599999986pt"/> is the height of the triangle treating <img alt="$\mathbf{e}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/2dfa893f7bf3c5c01ceb4788da782da0.svg" align="middle" width="13.31525249999999pt" height="14.611878600000017pt"/> as base, and
likewise for <img alt="$h_j$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6d22be1359e204374e6f0b45e318d561.svg" align="middle" width="15.57562379999999pt" height="22.831056599999986pt"/>. Rewriting the equation above, replacing one of the edge norms,
e.g.\ <img alt="$\|\mathbf{e}_i\|$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/d0cc021c5b0e6847b9929fd5c6ecd826.svg" align="middle" width="30.57556919999999pt" height="24.65753399999998pt"/>:
<p align="center"><img alt="$$&#10;-\frac{\|\mathbf{e}_j\|}{2A} \frac{\frac{h_j}{\sin{\alpha}_{ij}}}{2A} \cos {\alpha}_{ij}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/99e92f13a02d606f5ce1fd991794f828.svg" align="middle" width="148.3049766pt" height="42.926063400000004pt"/></p>


Combine the cosine and sine terms:
<p align="center"><img alt="$$&#10;-\frac{\|\mathbf{e}_j\|}{2A} \frac{h_j}{2A} \cot {\alpha}_{ij}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/bd842203a257378897dd2d2499081e97.svg" align="middle" width="127.30930079999999pt" height="34.7253258pt"/></p>


Finally, since <img alt="$\|\mathbf{e}_j\|h_j=2A$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/940aebca1f028d93419529cff94a740a.svg" align="middle" width="90.89233229999999pt" height="24.65753399999998pt"/>, our constant dot product of these
gradients in our triangle is:
<p align="center"><img alt="$$&#10;{\nabla} {\varphi}_i \cdot  {\nabla} {\varphi}_j = -\frac{\cot {\alpha}_{ij}}{2A}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/40cb87da7035dbcd178fcf1363089dfc.svg" align="middle" width="163.14031305pt" height="32.50746015pt"/></p>


Similarly, inside the other triangle <img alt="$T_{\beta}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e0563ab36735bded51644cf89825f1d2.svg" align="middle" width="17.65134359999999pt" height="22.465723500000017pt"/> incident
on nodes <img alt="$i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663225699999989pt" height="21.68300969999999pt"/> and <img alt="$j$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/36b5afebdba34564d884d347484ac0c7.svg" align="middle" width="7.710416999999989pt" height="21.68300969999999pt"/> with angle <img alt="${\beta}_{ij}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/45f257d25bb292b0f5b8a984f4efbe84.svg" align="middle" width="20.05337234999999pt" height="22.831056599999986pt"/> we have a constant dot
product:
<p align="center"><img alt="$$&#10;{\nabla} {\varphi}_i \cdot  {\nabla} {\varphi}_j = -\frac{\cot {\beta}_{ij}}{2B},&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/4f16a031aba6bcb3e69a09759c051f29.svg" align="middle" width="161.9226114pt" height="33.81208709999999pt"/></p>

where <img alt="$B$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/61e84f854bc6258d4108d08d4c4a0852.svg" align="middle" width="13.29340979999999pt" height="22.465723500000017pt"/> is the area <img alt="$T_{\beta}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e0563ab36735bded51644cf89825f1d2.svg" align="middle" width="17.65134359999999pt" height="22.465723500000017pt"/>.

Recall that <img alt="${\varphi}_i$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/965b9779f7d35773b987d1080ab428e9.svg" align="middle" width="15.404332349999988pt" height="14.15524440000002pt"/> and <img alt="${\varphi}_j$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/5e5d7941135a825fa8ab551e0eb76a4f.svg" align="middle" width="16.85794274999999pt" height="14.15524440000002pt"/> are only both nonzero inside these two
triangles, <img alt="$T_{\alpha}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/04fdc9ebe2a4f15ffdd732d0abd189e1.svg" align="middle" width="18.15209879999999pt" height="22.465723500000017pt"/> and <img alt="$T_{\beta}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e0563ab36735bded51644cf89825f1d2.svg" align="middle" width="17.65134359999999pt" height="22.465723500000017pt"/>.  So, since these constants are inside an
integral over area  we may write:
<p align="center"><img alt="$$&#10;\int\limits_\mathbf{S} {\nabla} {\varphi}_i \cdot  {\nabla} {\varphi}_j \;dA = &#10;\left.A{\nabla} {\varphi}_i \cdot  {\nabla} {\varphi}_j \right|_{T_{\alpha}} + \left.B{\nabla} {\varphi}_i \cdot  {\nabla} {\varphi}_j \right|_{T_{\beta}}&#10;=&#10;-\frac{1}{2} \left( \cot {\alpha}_{ij} + \cot {\beta}_{ij} \right).&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/9f8210836fb3ca3da88e84938171d18f.svg" align="middle" width="560.18018265pt" height="47.164758299999995pt"/></p>


## Mass matrix

Treated as an _operator_ (i.e., when used multiplied against a vector <img alt="$\mathbf{L}\mathbf{u}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/5435f08c809898be0b38d31c7232e2af.svg" align="middle" width="21.87202709999999pt" height="22.55708729999998pt"/>),
the Laplacian matrix <img alt="$\mathbf{L}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/80637df1ca7533740cc7b3fdd1ab540b.svg" align="middle" width="11.36979854999999pt" height="22.55708729999998pt"/> computes the local integral of the Laplacian of a
function <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/>. In the energy-based formulation of the smoothing problem this is
not an issue. If we used a similar FEM derivation for the _data term_ we would
get another sparse matrix <img alt="$\mathbf{M} \in  \mathbb{R}^{n \times  n}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e81a1a11e934fd54074805470b4f9ad1.svg" align="middle" width="76.43450429999999pt" height="26.17730939999998pt"/>:

<p align="center"><img alt="$$&#10;\int _\mathbf{S} (u-f)^{2} \;dA &#10;  = \int _\mathbf{S} {\sum}_{i=1}^n {\sum}_{j=1}^n {\varphi}_i\cdot {\varphi}_j (u_i-f_i) (u_j-f_j) \;dA =&#10;  (\mathbf{u}-\mathbf{f})^{\mathsf T} \mathbf{M} (\mathbf{u}-\mathbf{f}),&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e9b8736a9aa66b1921af8ad5a22571df.svg" align="middle" width="591.3941252999999pt" height="37.3519608pt"/></p>

where <img alt="$\mathbf{M}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg" align="middle" width="17.94511949999999pt" height="22.55708729999998pt"/> as an operator computes the local integral of a function's value
(i.e., <img alt="$\mathbf{M}\mathbf{u}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/5a6a7ddb47dbe4427df2767dd2254d3c.svg" align="middle" width="28.447348049999988pt" height="22.55708729999998pt"/>).

This matrix <img alt="$\mathbf{M}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg" align="middle" width="17.94511949999999pt" height="22.55708729999998pt"/> is often _diagonalized_ or _lumped_ into a diagonal matrix,
even in the context of FEM. So often we will simply set:

<p align="center"><img alt="$$&#10;M_{ij} = &#10;\begin{cases}&#10;  \frac13  {\sum}_{t=1}^m \begin{cases}&#10;  \text{Area}(t) &amp; \text{if triangle $t$ contains vertex $i$} \\&#10;  0 &amp; \text{otherwise}&#10;  \end{cases}&#10;  &amp; \text{if $i=j$}\\&#10;  0 &amp; \text{otherwise},&#10;\end{cases}&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/48f326d21a7a528a0f1aeb6a90755416.svg" align="middle" width="497.2228635pt" height="71.0142642pt"/></p>

for a mesh with <img alt="$m$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/0e51a2dede42189d77627c4d742822c3.svg" align="middle" width="14.433101099999991pt" height="14.15524440000002pt"/> triangles.

If we start directly with the continuous smoothing iteration equation, then we
have a point-wise equality. To fit in our integrated Laplacian <img alt="$\mathbf{L}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/80637df1ca7533740cc7b3fdd1ab540b.svg" align="middle" width="11.36979854999999pt" height="22.55708729999998pt"/> we should
convert it to a point-wise quantity. From a units perspective, we need to
divide by the local area. This would result in a discrete smoothing iteration
equation:

<p align="center"><img alt="$$&#10;\mathbf{u}^t = (\mathbf{I} - {\lambda}\mathbf{M}^{-1} \mathbf{L})\mathbf{u}^{t+1},&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/5d53591a59f665c794e0a2f859142b53.svg" align="middle" width="172.30560434999998pt" height="18.312383099999998pt"/></p>


where <img alt="$\mathbf{I} \in  \mathbb{R}^{n\times n}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/d4fecafea88cb244ca98be86cb2bb526.svg" align="middle" width="65.65830809999999pt" height="26.17730939999998pt"/> is the identity matrix. This equation is _correct_ but
the resulting matrix <img alt="$\mathbf{A} := \mathbf{I} - {\lambda}\mathbf{M}^{-1} \mathbf{L}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/4bdc180d4843f731cfac9a4b56ba3366.svg" align="middle" width="124.58859434999998pt" height="26.76175259999998pt"/> is not symmetric and thus slower
to solve against.

Instead, we could take the healthier view of requiring our smoothing iteration
equation to hold in a locally integrated sense. In this case, we replace mass
matrices on either side:

<p align="center"><img alt="$$&#10;\mathbf{M} \mathbf{u}^t = (\mathbf{M} - {\lambda}\mathbf{L})\mathbf{u}^{t+1}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/58bea700f9fd133dfde44969079927c0.svg" align="middle" width="165.4333428pt" height="18.312383099999998pt"/></p>


Now the system matrix <img alt="$\mathbf{A} := \mathbf{M} + {\lambda}\mathbf{L}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/c317e35fe20ee616ed6678d04169794b.svg" align="middle" width="99.77121329999999pt" height="22.831056599999986pt"/> will be symmetric and we can use
[Cholesky factorization](https://en.wikipedia.org/wiki/Cholesky_decomposition)
to solve with it.

### Laplace Operator is Intrinsic

The discrete Laplacian operator and its accompanying mass matrix are
_intrinsic_ operators in the sense that they _only_ depend on lengths. In
practical terms, this means we do not need to know _where_ vertices are
actually positioned in space (i.e., <img alt="$\mathbf{V}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/26eb59da31fb48cb17abfe4c6dc80375.svg" align="middle" width="14.554737449999989pt" height="22.55708729999998pt"/>). Rather we only need to know the
relative distances between neighboring vertices (i.e., edge lengths). We do not
even need to know which dimension this mesh is [living
in](https://en.wikipedia.org/wiki/Embedding).

This also means that applying a transformation to a shape that does not change
any lengths on the surface (e.g., bending a sheet of paper) will have no affect
on the Laplacian.

### Data denoising

For the data denoising application, our geometry of the domain is not changing
only the scalar function living upon it. We can build our discrete Laplacian
<img alt="$\mathbf{L}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/80637df1ca7533740cc7b3fdd1ab540b.svg" align="middle" width="11.36979854999999pt" height="22.55708729999998pt"/> and mass matrix <img alt="$\mathbf{M}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg" align="middle" width="17.94511949999999pt" height="22.55708729999998pt"/> and apply the above formula with a chosen <img alt="${\lambda}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/1b109d8b4484cf614f27126d788c510e.svg" align="middle" width="9.58908224999999pt" height="22.831056599999986pt"/>
parameter.

### Geometric smoothing

For geometric smoothing, the Laplacian operator (both <img alt="$\Delta $" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/545994479b8eb0934d80a97627e63f82.svg" align="middle" width="13.69867124999999pt" height="22.465723500000017pt"/> in the continuous
setting and <img alt="$\mathbf{L},\mathbf{M}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/a72b56a59ffa2be52208086ccb6f7562.svg" align="middle" width="36.62080289999999pt" height="22.55708729999998pt"/> in the discrete setting) depend on the geometry of the
surface <img alt="$\mathbf{S}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/>. So if the signal <img alt="$u$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/> is replaced with the positions of points on
the surface (say, <img alt="$\mathbf{V} \in  \mathbb{R}^{n\times 3}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/8fc2ed1e8ffd9008ef2661b89ca5b616.svg" align="middle" width="71.47064264999999pt" height="26.76175259999998pt"/> in the discrete case), then the smoothing
iteration update rule is a _non-linear_ function if we write it as:

<p align="center"><img alt="$$&#10;\mathbf{M}^{t+1} \mathbf{V}^t = (\mathbf{M}^{t+1} - {\lambda}\mathbf{L}^{t+1})\mathbf{V}^{t+1}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/c6f234133df83b61064c6cdda4c75eb8.svg" align="middle" width="240.8331519pt" height="18.312383099999998pt"/></p>


However, if we assume that small changes in <img alt="$\mathbf{V}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/26eb59da31fb48cb17abfe4c6dc80375.svg" align="middle" width="14.554737449999989pt" height="22.55708729999998pt"/> have a negligible effect on
<img alt="$\mathbf{L}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/80637df1ca7533740cc7b3fdd1ab540b.svg" align="middle" width="11.36979854999999pt" height="22.55708729999998pt"/> and <img alt="$\mathbf{M}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg" align="middle" width="17.94511949999999pt" height="22.55708729999998pt"/> then we can discretize _explicitly_ by computing <img alt="$\mathbf{L}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/80637df1ca7533740cc7b3fdd1ab540b.svg" align="middle" width="11.36979854999999pt" height="22.55708729999998pt"/> and <img alt="$\mathbf{M}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg" align="middle" width="17.94511949999999pt" height="22.55708729999998pt"/>
_before_ performing the update:

<p align="center"><img alt="$$&#10;\mathbf{M}^{t} \mathbf{V}^t = (\mathbf{M}^{t} - {\lambda}\mathbf{L}^{t}) \mathbf{V}^{t+1}.&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/8c168fd96735d211e2c39766d370606a.svg" align="middle" width="190.90141949999997pt" height="18.312383099999998pt"/></p>


### Why did my mesh disappear?

Repeated application of geometric smoothing may cause the mesh to "disappear".
Actually the updated vertex values are being set to
[NaNs](https://en.wikipedia.org/wiki/NaN) due to degenerate numerics. We are
rebuilding the discrete Laplacian at every new iteration, regardless of the
"quality" of the mesh's triangles. In particular, if a triangle tends to become
skinnier and skinnier during smoothing, what will happen to the cotangents of
its angles?

In "Can Mean-Curvature Flow Be Made Non-Singular?", Kazhdan et al. derive a new
type of geometric flow that is stable (so long as the mesh at time <img alt="$t=0$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/1c899e1c767eb4eac89facb5d1f2cb0d.svg" align="middle" width="36.07293689999999pt" height="21.18721440000001pt"/> is
reasonable). Their change is remarkably simple: do not update <img alt="$\mathbf{L}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/80637df1ca7533740cc7b3fdd1ab540b.svg" align="middle" width="11.36979854999999pt" height="22.55708729999998pt"/>, only update
<img alt="$\mathbf{M}$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/e6bb22a58889cb2e58f4fce2f3a80e02.svg" align="middle" width="17.94511949999999pt" height="22.55708729999998pt"/>.

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

<p align="center"><img alt="$$&#10;L_{ij} = \begin{cases}&#10;         \frac12  \cot{{\alpha}_{ij}} + \frac12  \cot{{\beta}_{ij}}  &amp; \text{if edge $ij$ exists} \\&#10;         - {\sum}_{j\ne i} L_{ij}                   &amp; \text{if $i = j$} \\&#10;         0                                &amp; \text{otherwise}&#10;         \end{cases}&#10;$$" src="https://cdn.jsdelivr.net/gh/alecjacobson/geometry-processing-smoothing@master/tex/b135eeed3e8770ac50d05128df513462.svg" align="middle" width="325.85501519999997pt" height="69.0417981pt"/></p>


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

