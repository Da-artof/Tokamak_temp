# Tokamak_temp
C program to calculate the temperature profile of a tokamak reactor.

The model equation to be solved is:

$$\frac{\partial T}{\partial t} = \frac{\partial}{\partial \theta}\left(Q<sup>11</sup>\frac{\partial T}{\partial \theta}\right) + \frac{\partial}{\partial \zeta}\left(Q<sup>22</sup>\frac{\partial T}{\partial \zeta}\right) + \frac{\partial}{\partial \theta}\left(Q<sup>12</sup>\frac{\partial T}{\partial \zeta}\right) + \frac{\partial}{\partial \zeta}\left(Q<sup>21</sup>\frac{\partial T}{\partial \theta}\right) - R\left(\theta,\zeta\right)T + S\left(\theta,\zeta\right)$$

For coefficients Q<sup>11</sup>, Q<sup>22</sup>, Q<sup>12</sup> and Q<sup>21</sup> that relate to the shape of the shell. R represents the connection of the shell to external heat sinks and S is the energy source due to the fusion plasma.

This is a 3D partial differential equation of order 2. The equation is solved using an implicit finite difference method, substituting in the finite difference expressions of the various derivatives, which gives us 5 point coupling. We reduce the 3D problem to a 2D problem, creating a matrix equation in space for each time step, solving this and using the computed values to find the values at the next time step.

For computational efficiency, we reorder indices to form banded matrices, which are less computationally expensive to invert. The outputted data is the computed value of T for each combination of $\theta$ and $\zeta$ at the final time step.
