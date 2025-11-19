# Model Predictive Control (MPC) — Theory (finite-horizon LQR)

This document provides a concise, self-contained theoretical description of a finite-horizon LQR-based Model Predictive Control (MPC). It contains the mathematical formulation, derivation of the optimal feedback law via dynamic programming (Riccati recursion), reference-tracking considerations, and practical numerical notes. No code or implementation references are included.

## 1. System model

Consider a discrete-time linear time-invariant dynamical system:

$$
x_{k+1} = A x_k + B u_k + w_k,
$$

where
- $x_k \in \mathbb{R}^n$ is the state at time step $k$,
- $u_k \in \mathbb{R}^m$ is the control input at time step $k$,
- $A \in \mathbb{R}^{n\times n}$ is the state transition matrix,
- $B \in \mathbb{R}^{n\times m}$ is the input matrix,
- $w_k$ is an optional disturbance (often ignored in the optimization stage).

In the optimal control problem below the disturbance $w_k$ is not explicitly optimized; it may be present in simulation or disturbance rejection analyses.

## 2. Finite-horizon quadratic cost

Define a finite prediction horizon $N > 0$ and the quadratic cost functional

$$
J(x_0, u_{0:N-1}) = \sum_{k=0}^{N-1} \left( x_k^\top Q x_k + u_k^\top R u_k \right) + x_N^\top Q_f x_N,
$$

where
- $Q \succeq 0$ is the state weighting matrix ($n\times n$),
- $R \succ 0$ is the control weighting matrix ($m\times m$),
- $Q_f \succeq 0$ is the terminal state weighting matrix.

The objective is to choose the control sequence $u_0,\dots,u_{N-1}$ minimizing $J$ subject to the linear dynamics.

## 3. Dynamic programming and quadratic value functions

Under quadratic cost and linear dynamics, the optimal cost-to-go (value function) at time $t$ can be written in quadratic form:

$$
V_t(x) = x^\top P_t x,
$$

for some symmetric positive semidefinite matrices $P_t \in \mathbb{R}^{n\times n}$, with terminal condition

$$
P_N = Q_f.
$$

The dynamic programming recursion for time $t$ (cost-to-go minimization for stage $t$) is

$$
V_t(x) = \min_u \left\{ x^\top Q x + u^\top R u + V_{t+1}(A x + B u) \right\}.
$$

Substituting the quadratic ansatz for $V_{t+1}$ yields:

$$
\begin{aligned}
J_t(x,u) &= x^\top Q x + u^\top R u + (A x + B u)^\top P_{t+1} (A x + B u) \\
&= x^\top (Q + A^\top P_{t+1} A) x + 2 u^\top (B^\top P_{t+1} A) x + u^\top (R + B^\top P_{t+1} B) u.
\end{aligned}
$$

Define intermediate matrices

$$
S_t := R + B^\top P_{t+1} B, \qquad G_t := B^\top P_{t+1} A,
$$

so the stage cost is

$$
J_t(x,u) = x^\top \tilde{P} x + 2 u^\top G_t x + u^\top S_t u,
$$

where $\tilde{P} = Q + A^\top P_{t+1} A$.

## 4. Optimal feedback and Riccati recursion

Minimizing the quadratic form in $u$ (for fixed $x$) gives the first-order condition:

$$
\frac{\partial J_t}{\partial u} = 2 S_t u + 2 G_t x = 0 \quad\Rightarrow\quad S_t u^\star = - G_t x.
$$

Assuming $S_t$ is invertible (or working with a suitable regularization), the optimal stage control is a linear state-feedback law:

$$
u_t^\star = - K_t x_t, \qquad K_t := S_t^{-1} G_t = (R + B^\top P_{t+1} B)^{-1} B^\top P_{t+1} A.
$$

Substituting $u_t^\star$ back into the expression for $J_t$ yields the Riccati update for $P_t$:

$$
P_t = Q + A^\top P_{t+1} A - A^\top P_{t+1} B K_t.
$$

Equivalently (explicitly showing the Schur complement form):

$$
P_t = Q + A^\top P_{t+1} A - A^\top P_{t+1} B (R + B^\top P_{t+1} B)^{-1} B^\top P_{t+1} A.
$$

The backward recursion proceeds from $t=N-1$ down to $t=0$, initializing with $P_N = Q_f$. The set \{K_0, K_1, \dots, K_{N-1}\} are the optimal time-varying feedback gains for the finite-horizon problem.

## 5. Receding-horizon implementation (MPC)

The receding-horizon implementation uses the finite-horizon solution in an online loop. At each sampling instant:

1. Measure or estimate the current state $x_t$.
2. Solve the finite-horizon optimal control problem (or use precomputed gains if dynamics are time-invariant) to obtain the sequence of optimal controls $u_t^\star, u_{t+1}^\star,\dots,u_{t+N-1}^\star$ (or their feedback gains). 
3. Apply only the first control in the sequence: $u_t = u_t^\star = -K_0 x_t$ (or the equivalent immediate feedback law obtained for time 0). 
4. At the next time step, obtain the new state and repeat the process (shift horizon forward).

This receding-horizon rule provides feedback and can handle model mismatch and disturbances by re-planning at each step.

## 6. Reference tracking and steady-state feedforward

To track a desired constant reference state $x_{\mathrm{ref}}$, it is common to decompose the control into a steady-state feedforward term $u_{\mathrm{ref}}$ plus feedback on the state error. A (steady-state) equilibrium condition for a constant reference is:

$$
x_{\mathrm{ref}} = A x_{\mathrm{ref}} + B u_{\mathrm{ref}} \quad\Longrightarrow\quad (I - A) x_{\mathrm{ref}} = B u_{\mathrm{ref}}.
$$

If this linear equation has a solution, one can pick $u_{\mathrm{ref}}$ satisfying it (for under/overdetermined cases, solve in the least-squares sense or use the pseudoinverse). Then define the error $\tilde{x}_t = x_t - x_{\mathrm{ref}}$. The regulator acts on the error:

$$
u_t = -K_0 \tilde{x}_t + u_{\mathrm{ref}}.
$$

This achieves setpoint regulation if $u_{\mathrm{ref}}$ is chosen so that $x_{\mathrm{ref}}$ is an equilibrium of the closed-loop dynamics with the feedforward term.

## 7. Numerical considerations and regularization

Several numerical issues commonly arise and are addressed in practical implementations:

- Inversion of $S_t = R + B^\top P_{t+1} B$: if $S_t$ is ill-conditioned or singular, it is typical to add a small regularization term $\epsilon I$ to ensure invertibility, or to solve linear systems via numerically stable factorizations (Cholesky, SVD).

- Symmetrization of Riccati matrices: floating-point arithmetic can introduce slight asymmetry. It is common to symmetrize $P_t$ by averaging with its transpose to maintain numerical symmetry.

- Solving for steady-state feedforward $u_{\mathrm{ref}}$: when $(I-A) x_{\mathrm{ref}} = B u_{\mathrm{ref}}$ is under- or over-determined, use regularized least-squares (solve $(B^\top B + \delta I) u = B^\top b$) or the Moore-Penrose pseudoinverse to find a numerically stable solution. Choose the regularization level $\delta$ small relative to the problem scale.

- Scaling of cost matrices $Q,R$: the relative magnitudes of $Q$ and $R$ determine the trade-off between state regulation and control effort. Poor scaling can produce ill-conditioned Riccati equations.

- Horizon length $N$: short horizons reduce computation and may degrade performance; long horizons approach infinite-horizon LQR but increase computation and can magnify model mismatch.

## 8. Special cases and relations

- Infinite-horizon LQR: As the horizon $N\to\infty$ for a stabilizable pair $(A,B)$ and with appropriate detectability conditions, the sequence $P_t$ (run backward) may converge to a fixed point $P$ satisfying the discrete algebraic Riccati equation (DARE). The corresponding constant gain $K$ is the infinite-horizon LQR gain.

- Unconstrained vs constrained MPC: The quadratic derivation above is for the unconstrained case (no input/state inequality constraints). Constrained MPC requires solving a Quadratic Program (QP) online; the unconstrained solution gives useful intuition and is used when constraints are absent or when fast, simple controllers are sufficient.

- Time-varying systems: If $A$ or $B$ vary with time, the Riccati recursion and gains $K_t$ become truly time-varying and must be recomputed for every horizon or handled by a time-varying LQR formulation.

## 9. Practical tips and robustification

- Add small positive regularization in $R$ (if necessary) to keep $S_t$ well-conditioned.
- Use numerically stable linear algebra (Cholesky decomposition for symmetric positive definite matrices, SVD for pseudoinverse, and solve linear systems instead of forming explicit inverses when possible).
- Monitor conditioning and the magnitudes of feedforward terms; avoid applying extremely large feedforward actions without sanity checks.
- When disturbances are significant, consider augmenting the state with disturbance estimates (observer) or using stochastic/LQG formulations.

## 10. Summary (compact list of central equations)

- Dynamics: $$x_{k+1} = A x_k + B u_k + w_k.$$ 
- Cost (finite horizon): $$J = \sum_{k=0}^{N-1} (x_k^\top Q x_k + u_k^\top R u_k) + x_N^\top Q_f x_N.$$ 
- Backward recursion (terminal): $$P_N = Q_f.$$ 
- Intermediate matrices: $$S_t = R + B^\top P_{t+1} B, \quad G_t = B^\top P_{t+1} A.$$ 
- Optimal feedback gain: $$K_t = S_t^{-1} G_t = (R + B^\top P_{t+1} B)^{-1} B^\top P_{t+1} A.$$ 
- Riccati update: $$P_t = Q + A^\top P_{t+1} A - A^\top P_{t+1} B K_t.$$ 
- Receding-horizon control law (applied online): $$u_t = -K_0 x_t\quad\text{(or)}\quad u_t = -K_0 (x_t - x_{\mathrm{ref}}) + u_{\mathrm{ref}}.$$ 

---

