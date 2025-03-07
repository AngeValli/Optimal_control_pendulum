from typing import Tuple
import parameters as p

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import scipy.optimize
from numdifftools import Jacobian, Hessian

def differential_system_minimise_time(z_0: np.ndarray[float]) -> np.ndarray[float]:
    """
    Compute the derivative of the state vector z
    - Params :
        z (np.array[float]) : State vector
    - Returns :
        z_dot (np.array[float]) : Derivative of the state vector
    """
    z_dot: np.ndarray[float] = np.zeros(4)
    z_dot[0] = z_0[1]
    u: float = -np.sign(z_0[3])*p.U_B
    z_dot[1] = -p.g/p.l * np.sin(z_0[0]) - p.k/(p.m*p.l**2) * z_0[1] + u/(p.m*p.l**2)
    z_dot[2] = z_0[3] * p.g/p.l * np.cos(z_0[0])
    z_dot[3] = -z_0[2] + p.k * z_0[3]/(p.m*p.l**2)
    return z_dot

def jacobian_criteria_to_minimise(params_: np.ndarray[float]) -> float:
    """
    Jacobian computation of the criteria
    """
    return Jacobian(lambda params_: criteria_minimise_time(params_)[0])(params_).ravel()

def hessian_criteria_to_minimise(params_: np.ndarray[float]) -> float:
    """
    Hessian computation of the criteria
    """
    return Hessian(lambda params_: criteria_minimise_time(params_)[0])(params_)

def criteria_minimise_time(params_: np.ndarray[float]) -> Tuple[float,np.ndarray[float], np.ndarray[float]]:
    """
    Compute the criterion and perform solving of the two boundaries value problem
    - Params :
        lambda_ (np.array[float]) : Initial co-states
    - Returns :
        error (float) : Error between the final state and the desired state
        z_opt (np.array[np.ndarray[float]]) : Optimal solution state values
        t_opt (np.array[float]) : Time steps corresponding to optimal solution
    """
    z_0: np.ndarray[float] = np.array([np.pi/2, 0, np.cos(params_[0]), np.sin(params_[0])])
    sol = solve_ivp(fun=lambda _, x: differential_system_minimise_time(x), t_span=[0, np.exp(params_[1])], y0=z_0, t_eval=np.arange(0, np.exp(params_[1]), p.TE))
    y_opt: np.ndarray[np.ndarray[float]] = sol.y
    t_opt: np.ndarray[float] = sol.t
    criteria: float = np.linalg.norm(y_opt[0:2, -1]) + np.linalg.norm(np.abs(y_opt[3, -1])-p.m*p.l**2/p.U_B)
    return criteria, y_opt, t_opt

def compute_energy_minimise_time(z_opt: np.ndarray[np.ndarray[float]], t_opt: np.ndarray[float]) -> Tuple[float, np.ndarray[float]]:
    """
    Compute the energy
    - Params :
        z_opt (np.array[np.ndarray[float]]) : Optimal solution state values
        t_opt (np.array[float]) : Time steps corresponding to optimal solution
    - Returns :
        energy (float) : energy, in Joules
        u_estim (np.array[float]) : Estimated control
    """
    u_estim: float = -np.sign(z_opt[3])*p.U_B
    energy: float = np.trapz(y=u_estim**2, x=t_opt)
    return energy, u_estim

def solve_minimise_time() -> Tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Solve the two boundaries value problem
    - Returns :
        z_opt (np.array[np.ndarray[float]]) : Optimal solution state values
        t_opt (np.array[float]) : Time steps corresponding to optimal solution
    """
    res = scipy.optimize.minimize(fun=lambda params_: criteria_minimise_time(params_)[0], x0=p.PARAMS_0, method="dogleg", jac=jacobian_criteria_to_minimise, hess=hessian_criteria_to_minimise, tol=p.TOL)
    _, z_opt, t_opt = criteria_minimise_time(params_=res.x)
    return z_opt, t_opt

def display_minimise_time(z_opt: np.ndarray[np.ndarray[float]], t_opt: np.ndarray[float]) -> None:
    """
    Display the results for the two boundaries value problem
    - Params :
        z_opt (np.array[np.ndarray[float]]) : Optimal solution state values
        t_opt (np.array[float]) : Time steps corresponding to optimal solution
    """
    energy: float
    u_estim: np.ndarray[float]
    energy, u_estim = compute_energy_minimise_time(z_opt, t_opt)

    fig1, ax1 = plt.subplots(1, 2, sharex=True)
    fig1.set_figheight(5)
    fig1.set_figwidth(15)
    ax1[0].plot(t_opt, z_opt[0,:], label="State $x_1(t)$")
    ax1[0].plot(t_opt, z_opt[1,:], label="State $x_2(t)$")
    ax1[0].set_title("States")
    ax1[0].set_xlabel("Time (s)")
    ax1[0].legend()
    ax1[1].plot(t_opt, z_opt[2,:], label="Costate $\\lambda_1(t)$")
    ax1[1].plot(t_opt, z_opt[3,:], label="Costate $\\lambda_2(t)$")
    ax1[1].set_title("Costates")
    ax1[1].set_xlabel("Time (s)")
    ax1[1].legend()
    plt.savefig("States_and_costates_minimise_time.png")

    fig2, ax2 = plt.subplots(1, 3)
    fig2.set_figheight(5)
    fig2.set_figwidth(15)
    ax2[0].plot(t_opt, z_opt[0,:], label="Angle $\\theta(t)$")
    ax2[0].plot(t_opt, z_opt[1,:], label="Speed $\\dot{\\theta}(t)$")
    ax2[0].set_title("Angle $\\theta(t)$ and speed $\\dot{\\theta}(t)$")
    ax2[0].set_xlabel("Time (s)")
    ax2[0].legend()
    ax2[1].plot(z_opt[0,:], z_opt[1,:])
    ax2[1].set_xlabel("Speed $\\dot{\\theta}(t)$")
    ax2[1].set_ylabel("Angle $\\theta(t)$")
    ax2[2].plot(t_opt, u_estim)
    ax2[2].set_title(f"Command $u(t)$  : Energy {energy:.2f} J")
    ax2[2].set_xlabel("Time (s)")
    plt.savefig("Command_minimise_time.png")

if __name__ == "__main__":
    z_opt: np.ndarray[np.ndarray[float]]
    t_opt: np.ndarray[float]
    z_opt, t_opt = solve_minimise_time()
    display_minimise_time(z_opt, t_opt)