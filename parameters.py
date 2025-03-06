from typing import List

import numpy as np

# Physical parameters
m: float = 1 # Mass
g: float = 10 # Gravity
l: float = 1 # Length
k: float = 0.1 # Coefficient of friction
# You can change it from 0.01 to 0.9 to observe different behaviors of the pendulum

# Solving the problem for minimum energy
T_MAX: int = 10 # Global final time
LAMBDA_0: np.ndarray[float] = np.array([0, 0]) # Initial costates
TOL: float = 1e-6 # Tolerance
TE: float = 1e-2 # Sampling period
T_SPAN: List[int] = [0, T_MAX] # Time span (to include T_MAX value in the solution)
T_EVAL: np.ndarray[float] = np.arange(0, T_MAX, TE) # Time evaluation
T_SPAN_VERIF: List[int] = [0, T_MAX] # Time span for verification (to exclude T_MAX value from the points to interpolate)
Z_0_VERIF: np.ndarray[float] = np.array([np.pi/2, 0]) # Initial state for verification

# Solving the problem for minimum time
U_B: int = 1 # Control bound
ALPHA_0: float = 0 # Initial costate value
T_MAX: int = 10
PARAMS_0: np.ndarray[float] = np.array([ALPHA_0, np.log(T_MAX)]) # Initial parameters