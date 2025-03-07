# Optimal_control_pendulum

Implementation of an optimal control strategy for a nonlinear system of a pendulum.

<img src="./pendulum.png" alt="pendulum" width="200"/>

The continuous-time nonlinear dynamic of the system can be described by the following ODE:

$$
m L^2 \ddot{\theta} = -m g L \sin(\theta)-k \dot{\theta} + u
$$

With $\theta$ the angular position of the pendulum, $u$ the command being the torque applied to the axis of rotation, $L$ the length of the rod, $g$ the force of gravity and $k$ a coefficient of fluid friction.

## Minimise the energy

The initial state is $(\theta, \dot{\theta}) = \big( \frac{\pi}{2}, 0 \big)$ and the goal is to reach the final state at $t=t_f$ fixed to the equilibrium of the pendulum $(\theta, \dot{\theta}) = (0, 0)$.

We want to minimise the energy used so the criteria is given by :

$$
\int_0^{t_f}u^2(t) dt
$$

Let define a state vector

$$
z =
\begin{bmatrix}
\theta \\
\dot{\theta}
\end{bmatrix} =
\begin{bmatrix}
z_1 \\
z_2
\end{bmatrix}
$$

The derivative has the following form :

$$
\frac{dz}{dt} = f(z)
$$

We obtain :

$$
\frac{dz}{dt} = \frac{d}{dt}
\begin{bmatrix}
z_1 \\
z_2
\end{bmatrix} =
\begin{bmatrix}
z_2 \\
-\frac{g}{L}\sin(z_1) - \frac{k}{mL^2}z_2 + \frac{u}{mL^2}
\end{bmatrix}
$$

The canonical form of the criteria is given by :

$$
\int_0^{t_f} q dt + \phi(z(t_f)) = \int_0^{t_f}u^2(t) dt
$$

By identification, we have $q=u^2$ and $\phi(z(t_f))=0$. From the transversality conditions we obtain :

$$
\psi(t_f, z(t_f)) =
\begin{bmatrix}
\psi_1 \\
\psi_2
\end{bmatrix} =
\begin{bmatrix}
z_1(t_f) \\
z_2(t_f)
\end{bmatrix} =
\begin{bmatrix}
0 \\
0
\end{bmatrix}
$$

Let 

$$
\lambda =
\begin{bmatrix}
\lambda_1 \\
\lambda_2
\end{bmatrix}
$$

With $\lambda_1$ the costate associated to $x_1$ and $\lambda_2$ the costate associated to $x_2$, we define the hamiltonian $H$ such as :

$$
H = q + \lambda^T f(z) = u^2 + \lambda_1 z_2 + \lambda_2 \bigg( -\frac{g}{L}\sin(z_1) - \frac{k}{mL^2}z_2 + \frac{u}{mL^2} \bigg)
$$

We write the Euler-Lagrange equations of the hamiltonian. First for the costates :

$$
\dot{\lambda} = - \frac{\partial H}{\partial z} \Longleftrightarrow
\begin{cases} \dot{\lambda}_1 = - \frac{\partial H}{\partial z_1} = \frac{\lambda_2 g}{L} \cos(z_1) \\
\dot{\lambda}_2 = - \frac{\partial H}{\partial z_2} = -\lambda_1 + \frac{k \lambda_2}{mL^2}
\end{cases}
$$

Then for the command $u(t)$, the optimal command $u^*$ is given by the condition:

$$
\frac{\partial H(u^\*)}{\partial u} = 2u^* + \frac{\lambda_2}{mL^2} = 0 \Longleftrightarrow u^* = - \frac{\lambda_2}{2mL^2}
$$

The transversality conditions gives the Lagrange multipliers

$$
\nu =
\begin{bmatrix}
\nu_1 \\
\nu_2
\end{bmatrix}
$$

such as

$$
\lambda(t_f) = \frac{\partial}{\partial z(t_f)}(\phi + \nu^T\psi) = 0 \Longleftrightarrow
\begin{cases}
\lambda_1(t_f) = \nu_1 \\
\lambda_2(t_f) = \nu_2
\end{cases}
$$

Finally, we obtain the following 2 boundaries value problem to find the initial values of the costates $\lambda_1(0)$ and $\lambda_2(0)$.

$$
\begin{cases}
\dot{z}_1 = z_2 \\
\dot{z}_2 = - \frac{g}{L}\sin(z_1)- \frac{k}{mL^2}z_2+\frac{u^*}{mL^2} \\
\dot{\lambda}_1 = \frac{\lambda_2 g}{L}\cos(z_1) \\
\dot{\lambda}_2 = - \lambda_1 + \frac{k \lambda_2}{mL^2} \\
x_1(0) = \frac{\pi}{2}, \quad x_1(t_f) = 0, \quad x_2(0) = 0, \quad x_2(t_f) = 0
\end{cases}
$$

From Euler-Lagrange equation on the optimal command $u^*(t)$, the system can be written :

$$
\begin{cases}
\dot{z}_1 = z_2 \\
\dot{z}_2 = - \frac{g}{L}\sin(z_1)- \frac{k}{mL^2}z_2-\frac{\lambda_2}{2m^2L^4} \\
\dot{\lambda}_1 = \frac{\lambda_2 g}{L}\cos(z_1) \\
\dot{\lambda}_2 = - \lambda_1 + \frac{k \lambda_2}{mL^2} \\
x_1(0) = \frac{\pi}{2}, \quad x_1(t_f) = 0, \quad x_2(0) = 0, \quad x_2(t_f) = 0
\end{cases}
$$

The final values of the costates $\lambda_1(t_f)$ and $\lambda_2(t_f)$ corresponds to Lagrange multipliers. The unknown parameters are the initial values of the costates. We want to find $\lambda_1(t_f)$ and $\lambda_2(t_f)$ such as we obtain $x_1(t_f)=x_2(t_f)=0$. From the shooting method, we generate multiples values for

$$
\lambda(0) =
\begin{bmatrix}
\lambda_1(0) \\
\lambda_2(0)
\end{bmatrix}
$$

such that

$$
\min_{\lambda(0)}
\begin{pmatrix}
\begin{Vmatrix}
z_1(t_f) \\
z_2(t_f)
\end{Vmatrix}
\end{pmatrix}
$$

and use a Quasi-Newton approach to solve this problem.

## Minimise the time

We suppose to have a bounded command $\forall t, \lvert u(t) \rvert \leq u_b$. In this problem, the time $t_f$ is not fixed as it is for the previous problem, and we the same state vector $z$ and state equation $\frac{dz}{dt} = f(z)$.

The criteria can be written as :

$$
\min(t_f) = \min \bigg( \int_0^{t_f} dt \bigg)
$$

The identification with the canonical form :

$$
\int_0^{t_f} q dt + \phi(z(t_f)) = \int_0^{t_f} dt
$$

so $q=1$ and $\phi(z(t_f)) = 0$. Same as before, the transversality conditions give :

$$
\psi(t_f, z(t_f)) =
\begin{bmatrix}
\psi_1 \\
\psi_2
\end{bmatrix} =
\begin{bmatrix}
z_1(t_f) \\
z_2(t_f)
\end{bmatrix} =
\begin{bmatrix}
0 \\
0
\end{bmatrix}
$$

And the hamiltonian is :

$$
H = q + \lambda^T f(z) = 1 + \lambda_1 z_2 + \lambda_2 \bigg( -\frac{g}{L}\sin(z_1) - \frac{k}{mL^2}z_2 + \frac{u}{mL^2} \bigg)
$$

So the Euler-Lagrange equations for the costates remains the same.

$$
\dot{\lambda} = - \frac{\partial H}{\partial z} \Longleftrightarrow
\begin{cases} \dot{\lambda}_1 = - \frac{\partial H}{\partial z_1} = \frac{\lambda_2 g}{L} \cos(z_1) \\
\dot{\lambda}_2 = - \frac{\partial H}{\partial z_2} = -\lambda_1 + \frac{k \lambda_2}{mL^2}
\end{cases}
$$

But the condition on the command $u$ becomes :

$$
\frac{\partial H(u^*)}{\partial u} = \frac{\lambda_2}{mL^2} = 0 \neq f(u)
$$

By hypothesis, the command $u$ is bounded. Therefore, from the definition of the hamiltonian, we remark it is linear in $u$. As we also have a saturated command, we can conclude the optimal command $u^*$ is bang-bang such as :

$$
\begin{align*}
& \lambda_2 > 0 \implies u^* = -u_b \\
& \lambda_2 < 0 \implies u^* = u_b
\end{align*}
$$

The goal is to obtain the pendulum to its equilibrium position $(\theta, \dot{\theta}) = (0, 0)$ so $x_1(t_f) = x_2(t_f)=0$

From the transversality conditions, we are in the case where the final time is free so :

$$
\lambda(t_f)=\frac{\partial}{\partial(z_{t_f})}(\phi+\nu^T\psi) = 0 \implies
\begin{cases}
\lambda_1(t_f) = \nu_1 \\
\lambda_2(t_f) = \nu_2
\end{cases}
$$

Then, at the final time we have :

$$
\begin{align*}
&H(x(t_f), u(t_f), \lambda(t_f), t_f) + \frac{\partial}{\partial t_f}\bigg( \phi(t_f, x(t_f))+\psi(t_f, x(t_f))^T\nu \bigg) = 0 \\
&\implies H(x(t_f), u(t_f), \lambda(t_f), t_f) = 1 + \frac{u^* (t_f)}{mL^2} \lambda_2(t_f) = 0 \\
\end{align*}
$$

Finally :

$$
1-\frac{\text{sign}(\lambda_2(t_f))u_b}{mL^2}\lambda_2(t_f) = 0 \Longleftrightarrow \lvert \lambda_2(t_f) \rvert =\frac{mL^2}{u_b}
$$

The 2 boundaries value problem can be written as :

$$
\begin{cases}
\dot{z}_1 = z_2 \\
\dot{z}_2 = - \frac{g}{L}\sin(z_1)- \frac{k}{mL^2}z_2+\frac{u^*}{mL^2} \\
\dot{\lambda}_1 = \frac{\lambda_2 g}{L}\cos(z_1) \\
\dot{\lambda}_2 = - \lambda_1 + \frac{k \lambda_2}{mL^2} \\
x_1(0) = \frac{\pi}{2}, \quad x_1(t_f) = 0, \quad x_2(0) = 0, \quad x_2(t_f) = 0 \\
\lvert \lambda_2(t_f) \rvert = \frac{mL^2}{u_b}
\end{cases}
$$

From Euler-Lagrange equation on the optimal command $u^*(t)$, the system can be written :

$$
\begin{cases}
\dot{z}_1 = z_2 \\
\dot{z}_2 = - \frac{g}{L}\sin(z_1)- \frac{k}{mL^2}z_2+ \text{sign}(\lambda_2)\frac{u_b}{2m^2L^4} \\
\dot{\lambda}_1 = \frac{\lambda_2 g}{L}\cos(z_1) \\
\dot{\lambda}_2 = - \lambda_1 + \frac{k \lambda_2}{mL^2} \\
x_1(0) = \frac{\pi}{2}, \quad x_1(t_f) = 0, \quad x_2(0) = 0, \quad x_2(t_f) = 0 \\
\lvert \lambda_2(t_f) \rvert = \frac{mL^2}{u_b}
\end{cases}
$$

As in the previous problem, the final values of the costates $\lambda_1(t_f)$ and $\lambda_2(t_f)$ corresponds to Lagrange multipliers. The unknown parameters are the initial values of the costates. We want to find $\lambda_1(t_f)$ and $\lambda_2(t_f)$ such as we obtain $x_1(t_f)=x_2(t_f)=0$ and $H(t_f)=0$. From the shooting method, we generate multiples values for

$$
\lambda(0) =
\begin{bmatrix}
\lambda_1(0) \\
\lambda_2(0)
\end{bmatrix}
$$

such that

$$
\min_{t_f,\lambda(0)}
\begin{pmatrix}
\begin{Vmatrix}
z_1(t_f) \\
z_2(t_f) \\
H(t_f)
\end{Vmatrix}
\end{pmatrix}
$$

To simplify the expression, we can rewrite the costates equations such that :

$$
\begin{cases}
\dot{\lambda}_1 = \frac{\lambda_2 g}{L}\cos(z_1) \\
\dot{\lambda}_2 = - \lambda_1 + \frac{k \lambda_2}{mL^2}
\end{cases} \Longleftrightarrow
\dot{\lambda} = A(t) \lambda
$$

with

$$
A(t) =
\begin{pmatrix}
0 & \frac{g}{L}\cos(z_1) \\
-1 & \frac{k}{mL^2}
\end{pmatrix}
$$

Therefore any solution to this equation multiplied by a constant with respect to the time is a solution. We can look for a normalised solution on the unit circle

$$
\lVert \lambda(0) \rVert = 1
$$

We can represent this unit vector by its polar coordinates, meaning there exists an angle $\alpha$ such as :

$$
\begin{cases}
\lambda_1(0) = \cos(\alpha) \\
\lambda_2(0) = \sin(\alpha)
\end{cases}
$$

So we look for the minimisation problem over the angle $\alpha$ and use a Dogleg approach.
