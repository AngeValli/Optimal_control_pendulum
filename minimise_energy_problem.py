from typing import Tuple
import parameters as p

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import scipy.optimize

def differential_system_minimise_energy(z_0: np.ndarray[float]) -> np.ndarray[float]:
    """
    Compute the derivative of the state vector z
    - Params :
        z (np.array[float]) : State vector
    - Returns :
        z_dot (np.array[float]) : Derivative of the state vector
    """
    z_dot: np.ndarray[float] = np.zeros(4)
    z_dot[0] = z_0[1]
    u: float = -z_0[3]/(2*p.m*p.l**2)
    z_dot[1] = -p.g/p.l * np.sin(z_0[0]) - p.k/(p.m*p.l**2) * z_0[1] + u/(p.m*p.l**2)
    z_dot[2] = z_0[3] * p.g/p.l * np.cos(z_0[0])
    z_dot[3] = -z_0[2] + p.k * z_0[3]/(p.m*p.l**2)
    return z_dot

def criteria_minimise_energy(lambda_: np.ndarray[float]) -> Tuple[float,np.ndarray[float], np.ndarray[float]]:
    """
    Compute the criterion and perform solving of the two boundaries value problem
    - Params :
        lambda_ (np.array[float]) : Initial co-states
    - Returns :
        error (float) : Error between the final state and the desired state
        z_opt (np.array[np.ndarray[float]]) : Optimal solution state values
        t_opt (np.array[float]) : Time steps corresponding to optimal solution
    """
    z_0: np.ndarray[float] = np.array([np.pi/2, 0, lambda_[0], lambda_[1]])
    sol = solve_ivp(fun=lambda _, x: differential_system_minimise_energy(x), t_span=p.T_SPAN, y0=z_0, t_eval=p.T_EVAL)
    y_opt: np.ndarray[np.ndarray[float]] = sol.y
    t_opt: np.ndarray[float] = sol.t
    error: float = np.linalg.norm(y_opt[0:2, -1])
    return error, y_opt, t_opt

def solve_minimise_energy() -> Tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Solve the two boundaries value problem
    - Returns :
        z_opt (np.array[np.ndarray[float]]) : Optimal solution state values
        t_opt (np.array[float]) : Time steps corresponding to optimal solution
    """
    lambda_opt: np.ndarray[float]
    lambda_opt = scipy.optimize.fmin_bfgs(f=lambda lambda_: criteria_minimise_energy(lambda_)[0], x0=p.LAMBDA_0, gtol=p.TOL, xrtol=p.TOL)
    _, z_opt, t_opt = criteria_minimise_energy(lambda_=lambda_opt)
    return z_opt, t_opt

def compute_energy_minimise_energy(z_opt: np.ndarray[np.ndarray[float]], t_opt: np.ndarray[float]) -> Tuple[float, np.ndarray[float]]:
    """
    Compute the energy
    - Params :
        z_opt (np.array[np.ndarray[float]]) : Optimal solution state values
        t_opt (np.array[float]) : Time steps corresponding to optimal solution
    - Returns :
        energy (float) : energy, in Joules
        u_estim (np.array[float]) : Estimated control
    """
    u_estim: float = -z_opt[3, :]/(2*p.m*p.l**2)
    energy: float = np.trapz(y=u_estim**2, x=t_opt)
    return energy, u_estim

def display_minimise_energy(z_opt: np.ndarray[np.ndarray[float]], t_opt: np.ndarray[float]) -> None:
    """
    Display the results for the two boundaries value problem
    - Params :
        z_opt (np.array[np.ndarray[float]]) : Optimal solution state values
        t_opt (np.array[float]) : Time steps corresponding to optimal solution
    """
    energy: float
    u_estim: np.ndarray[float]
    energy, u_estim = compute_energy_minimise_energy(z_opt, t_opt)

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
    plt.savefig("States_and_costates_minimise_energy.png")

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
    plt.savefig("Command_minimise_energy.png")

if __name__ == "__main__":
    z_opt: np.ndarray[np.ndarray[float]]
    t_opt: np.ndarray[float]
    z_opt, t_opt = solve_minimise_energy()
    display_minimise_energy(z_opt, t_opt)