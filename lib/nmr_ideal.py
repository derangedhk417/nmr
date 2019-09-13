# This file contains definitions of expressions related to Nucelar Magnetic
# Resonance. It also includes definitions of related things, such as Pauli Spin
# Matrices and eigenspinors, becuase they are useful for calculating 
# expecation values. 
#
# The definitions in this file pertain to the idealized case where the 
# Time-Dependent Schrodinger equation can be solved perfectly. This file does
# not include calculations that account for transient effects.

from sympy                 import *
from IPython.display       import display, Math, Latex
from sympy.physics.quantum import Operator, Ket, Bra

# =============================================================================
# Principal Definitions
# =============================================================================

# Symbols
omega_0, omega, Omega  = symbols(r"\omega_0 \omega \Omega", real=True, positive=True)
omega_p                = symbols(r"\omega^{\prime}",        real=True, positive=True)
B_0, B_rf              = symbols(r"B_0 B_{rf}",             real=True, positive=True)
a_0, b_0               = symbols(r"a_0 b_0")
a, b                   = symbols(r"a b")
B, ih, jh, kh          = symbols(r"B \hat{i} \hat{j} \hat{k}")
theta                  = symbols(r"\theta", real=True, positive=True)
a_dot, b_dot           = symbols(r"\dot{a} \dot{b}")
t                      = symbols(r"t", real=True, positive=True)
gamma                  = symbols(r"\gamma", real=True, positive=True)

# Definitions Made for Clarity of Expression
expr_omega_p = sqrt((omega - omega_0)**2 + Omega**2)
expr_omega_0 = gamma*B_0
expr_Omega  =  gamma*B_rf
expr_a_0     = sin(theta / 2)
expr_b_0     = cos(theta / 2)

# The magnetic field
expr_B = B_rf*cos(omega*t)*ih + -B_rf*sin(omega*t)*jh + B_0*kh

# Differential Equation System
expr_a_dot = (I/2)*(Omega*exp(I*omega*t)*b + omega_0*a)
expr_b_dot = (I/2)*(Omega*exp(-I*omega*t)*a - omega_0*b)

# Solution to the System
expr_a = (a_0*cos(omega_p*t/2) + (I/omega_p)*(a_0*(omega_0 - omega) + b_0*Omega)*sin(omega_p*t/2))*exp(I*omega*t/2)
expr_b = (b_0*cos(omega_p*t/2) + (I/omega_p)*(b_0*(omega - omega_0) + a_0*Omega)*sin(omega_p*t/2))*exp(-I*omega*t/2)

# Complete Expressions, Without use of Definitions
complete_a = expr_a.subs(a_0, expr_a_0).subs(b_0, expr_b_0)
complete_a = complete_a.subs(omega_p, expr_omega_p)
complete_a = complete_a.subs(omega_0, expr_omega_0)
complete_a = complete_a.subs(Omega, expr_Omega)

complete_b = expr_b.subs(a_0, expr_a_0).subs(b_0, expr_b_0)
complete_b = complete_b.subs(omega_p, expr_omega_p)
complete_b = complete_b.subs(omega_0, expr_omega_0)
complete_b = complete_b.subs(Omega, expr_Omega)

# =============================================================================
# Spin Related Definitions
# =============================================================================

hbar = symbols(r"\hbar", real=True, positive=True)
sigma_x  = Operator(r"\mathbf{\sigma}_x")
sigma_y  = Operator(r"\mathbf{\sigma}_y")
sigma_z  = Operator(r"\mathbf{\sigma}_z")

expr_sigma_x = Matrix([[0, 1],  [1, 0]])
expr_sigma_y = Matrix([[0, -I], [I, 0]])
expr_sigma_z = Matrix([[1, 0],  [0, -1]])

S_x = Operator("S_x")
S_y = Operator("S_y")
S_z = Operator("S_z")

expr_S_x = (hbar/2)*expr_sigma_x
expr_S_y = (hbar/2)*expr_sigma_y
expr_S_z = (hbar/2)*expr_sigma_z

# =============================================================================
# Integrity Checks
# =============================================================================

# This code will check the above expressions to make sure that they solve the
# simplified system of equations, the full Time-Dependent Schrodinger equation
# and are orthonormal.

expr_ortho      = a*conjugate(a) + b*conjugate(b)
expr_ortho_eval = expr_ortho.subs(a, complete_a).subs(b, complete_b).simplify()

left_a  = Derivative(a, t)
right_a = expr_a_dot

left_b  = Derivative(b, t)
right_b = expr_b_dot

left_a_done  = left_a.subs(a, complete_a).doit()
right_a_done = right_a.subs(a, complete_a).subs(b, complete_b)
right_a_done = right_a_done.subs(Omega, expr_Omega)
right_a_done = right_a_done.subs(omega_0, expr_omega_0)
right_a_done = right_a_done.subs(a_0, expr_a_0).subs(b_0, expr_b_0)

left_b_done  = left_b.subs(b, complete_b).doit()
right_b_done = right_b.subs(a, complete_a).subs(b, complete_b)
right_b_done = right_b_done.subs(Omega, expr_Omega)
right_b_done = right_b_done.subs(omega_0, expr_omega_0)
right_b_done = right_b_done.subs(a_0, expr_a_0).subs(b_0, expr_b_0)

H = Operator(r"H")
expr_full_solution = Matrix([[expr_a], [expr_b]])
expr_sym_solution  = Matrix([[a], [b]])
expr_hamiltonian   = -(gamma * (hbar / 2))*(B_rf*(cos(omega*t)*expr_sigma_x - sin(omega*t)*expr_sigma_y) + B_0*expr_sigma_z)

expr_disp_left_schrodinger  = I*hbar*Derivative((UnevaluatedExpr(expr_sym_solution)), t)
expr_disp_right_schrodinger = UnevaluatedExpr(expr_hamiltonian)*UnevaluatedExpr(expr_sym_solution)

expr_left_a_simple = I*hbar*Derivative(a, t)
expr_left_b_simple = I*hbar*Derivative(b, t)

expr_right_a_simple = expr_disp_right_schrodinger.doit()[0]
expr_right_b_simple = expr_disp_right_schrodinger.doit()[1]

# =============================================================================
# Under Resonance
# =============================================================================

# Rewrite the equations for the resonance condition, but keep them compact.
expr_a_res = complete_a.subs(omega, B_0*gamma).subs(B_rf*gamma, Omega).subs(B_0*gamma, omega_0)
expr_b_res = complete_b.subs(omega, B_0*gamma).subs(B_rf*gamma, Omega).subs(B_0*gamma, omega_0)


# =============================================================================
# Expectation Values Under Resonance
# =============================================================================

S_x_exp    = symbols(r"\left<S_x\right>")
S_y_exp    = symbols(r"\left<S_y\right>")
S_z_exp    = symbols(r"\left<S_z\right>")
spin_state = symbols(r"\chi\left(t\right)")

expr_spin_generic = Matrix([[a], [b]])
expr_spin_state   = Matrix([[expr_a_res], [expr_b_res]])

# These are for Display Purposes
expr_S_x_exp_s     = UnevaluatedExpr(transpose(conjugate(expr_spin_generic)))
expr_S_x_exp_s    *= S_x*UnevaluatedExpr(expr_spin_generic)

expr_S_y_exp_s    = UnevaluatedExpr(transpose(conjugate(expr_spin_generic)))
expr_S_y_exp_s   *= S_y*UnevaluatedExpr(expr_spin_generic)

expr_S_z_exp_s    = UnevaluatedExpr(transpose(conjugate(expr_spin_generic)))
expr_S_z_exp_s   *= S_z*UnevaluatedExpr(expr_spin_generic)

# These will be evaluated
expr_S_x_exp     = trigsimp(((transpose(conjugate(expr_spin_state)))*expr_S_x*(expr_spin_state))[0].simplify())
expr_S_y_exp     = trigsimp(((transpose(conjugate(expr_spin_state)))*expr_S_y*(expr_spin_state))[0].simplify())
expr_S_z_exp     = trigsimp(((transpose(conjugate(expr_spin_state)))*expr_S_z*(expr_spin_state))[0].simplify())

# =============================================================================
# Thermal Equilibrium
# =============================================================================

N_down = symbols(r"N_down")
N_up   = symbols(r"N_up")
k_b    = symbols(r"k_b", real=True)
T      = symbols(r"T", real=True, positive=True)

expr_ratio        = exp(-(omega_0*hbar)/(k_b*T))
equil_solutions   = solve((expr_b_0 / expr_a_0) - expr_ratio, theta)
theta_equilibrium = equil_solutions[0]
expr_lim          = Limit(theta_equilibrium, T, oo)


# =============================================================================
# Expectation Values In Thermal Equilibrium
# =============================================================================

epsilon = symbols(r"\epsilon", real=True, positive=True)
theta_equil_approx = (pi / 2) + epsilon

expr_epsilon = (pi / 2) - theta_equilibrium

expr_a_0_equil = expr_a_0.subs(theta, theta_equil_approx)
expr_b_0_equil = expr_b_0.subs(theta, theta_equil_approx)

expr_S_x_exp_equil = expr_S_x_exp.subs(theta, theta_equil_approx)
expr_S_y_exp_equil = expr_S_y_exp.subs(theta, theta_equil_approx)
expr_S_z_exp_equil = expr_S_z_exp.subs(theta, theta_equil_approx)

expr_S_x_exp_equil_approx = expr_S_x_exp_equil.subs(sin(epsilon), epsilon)
expr_S_x_exp_equil_approx = expr_S_x_exp_equil_approx.subs(cos(epsilon), 1)

expr_S_y_exp_equil_approx = expr_S_y_exp_equil.subs(sin(epsilon), epsilon)
expr_S_y_exp_equil_approx = expr_S_y_exp_equil_approx.subs(cos(epsilon), 1)

expr_S_z_exp_equil_approx = expr_S_z_exp_equil.subs(sin(epsilon), epsilon)
expr_S_z_exp_equil_approx = expr_S_z_exp_equil_approx.subs(cos(epsilon), 1)

# This function will return nothing under normal circumstances. If a fundamental condition
# necessary for the expressions above to be correct is not satisfied, this will throw an
# exception.
def check_integrity():
	expr_left_a = simplify(expr_left_a_simple.subs(a, complete_a).subs(b, complete_b).doit())
	expr_left_b = simplify(expr_left_b_simple.subs(a, complete_a).subs(b, complete_b).doit())

	expr_right_a = simplify(expr_right_a_simple.subs(a, complete_a).subs(b, complete_b))
	expr_right_b = simplify(expr_right_b_simple.subs(a, complete_a).subs(b, complete_b))

	solve_a = simplify(expr_left_a - expr_right_a)
	solve_b = simplify(expr_left_b - expr_right_b)

	if expr_ortho_eval != 1:
		msg = "Othornormality is broken."
		raise Exception(msg)

	if simplify(left_a_done - right_a_done) != 0:
		msg = "The spin up part of the differential equation is not satisfied."
		raise Exception(msg)

	if simplify(left_b_done - right_b_done) != 0:
		msg = "The spin down part of the differential equation is not satisfied."
		raise Exception(msg)

	if solve_a != 0:
		msg = "The spin up part of the Schrodinger equation is not satisfied."
		raise Exception(msg)

	if solve_b != 0:
		msg = "The spin down part of the Schrodinger equation is not satisfied."
		raise Exception(msg)

















