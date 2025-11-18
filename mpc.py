"""mpc.py

Simple Model Predictive Controller implementation using finite-horizon LQR
(Riccati recursion). This implementation does not solve constrained QPs; it
computes the optimal unconstrained sequence for the quadratic cost
and applies the first control (receding horizon). Useful for controlling
linear networks x_{k+1} = A x_k + B u_k.

API:
  MPCController(A, B, Q, R, N, Qf=None)
    - compute finite-horizon LQR gains for horizon N and expose a .control(x)
      method that returns the receding-horizon control u = -K_0 x.

Functions:
  finite_horizon_lqr(A, B, Q, R, N, Qf=None)

Example usage is provided in the __main__ block. If A/B cannot be loaded
from workspace files, a small random example is used.
"""

from typing import List, Tuple, Optional
import numpy as np
import os
import csv


def finite_horizon_lqr(A: np.ndarray,
                       B: np.ndarray,
                       Q: np.ndarray,
                       R: np.ndarray,
                       N: int,
                       Qf: Optional[np.ndarray] = None) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Compute finite-horizon LQR gains using backward Riccati recursion.

    Solves the discrete-time finite-horizon LQR problem for the linear
    time-invariant system x_{k+1}=A x_k + B u_k with cost
      sum_{k=0}^{N-1} (x_k' Q x_k + u_k' R u_k) + x_N' Qf x_N

    Returns (Ks, Ps) where Ks[t] is the state-feedback gain at time t
    (u_t = -Ks[t] x_t) and Ps[t] is the Riccati matrix at time t.

    Notes:
      - No constraints are supported.
      - Matrices Q and R must be symmetric positive semidefinite/definite.
    """
    n = A.shape[0]
    m = B.shape[1]

    if Qf is None:
        Qf = Q

    # Pre-allocate lists
    Ps: List[np.ndarray] = [np.zeros((n, n)) for _ in range(N + 1)]
    Ks: List[np.ndarray] = [np.zeros((m, n)) for _ in range(N)]

    Ps[N] = Qf.copy()

    # Backward Riccati recursion
    for t in range(N - 1, -1, -1):
        Pnext = Ps[t + 1]
        # Compute optimal gain K_t = (R + B' Pnext B)^{-1} B' Pnext A
        S = R + B.T @ Pnext @ B
        # Ensure S is symmetric
        S = 0.5 * (S + S.T)
        try:
            Sinv = np.linalg.inv(S)
        except np.linalg.LinAlgError:
            # Regularize small eigenvalues if needed
            eps = 1e-8
            Sinv = np.linalg.inv(S + eps * np.eye(m))

        Kt = Sinv @ (B.T @ Pnext @ A)
        Ks[t] = Kt

        # Riccati update
        Ps[t] = Q + A.T @ Pnext @ A - A.T @ Pnext @ B @ Kt
        # Symmetrize for numerical stability
        Ps[t] = 0.5 * (Ps[t] + Ps[t].T)

    return Ks, Ps


class MPCController:
    """Receding-horizon controller using finite-horizon LQR gains.

    Usage:
      mpc = MPCController(A, B, Q, R, N)
      u = mpc.control(x)

    The controller precomputes a sequence of time-indexed gains for the
    horizon N and returns the first control. For time-varying A/B
    one can re-create the controller each step or extend this class.
    """

    def __init__(self,
                 A: np.ndarray,
                 B: np.ndarray,
                 Q: np.ndarray,
                 R: np.ndarray,
                 N: int = 10,
                 Qf: Optional[np.ndarray] = None):
        self.A = A
        self.B = B
        self.Q = Q
        self.R = R
        self.N = int(N)
        self.Qf = Qf if Qf is not None else Q

        # Validate shapes
        n = A.shape[0]
        assert A.shape[1] == n, "A must be square"
        assert B.shape[0] == n, "B must have same number of rows as A"
        assert Q.shape == (n, n), "Q must be n x n"
        assert R.shape[0] == R.shape[1], "R must be square"

        # Precompute gains
        self.Ks, self.Ps = finite_horizon_lqr(A, B, Q, R, self.N, self.Qf)

    def control(self, x: np.ndarray) -> np.ndarray:
        """Compute control action for current state x using first feedback gain.

        Returns u = -K_0 x
        """
        return self.control_with_ref(x)

    def control_with_ref(self, x: np.ndarray, x_ref: Optional[np.ndarray] = None, u_ref: Optional[np.ndarray] = None) -> np.ndarray:
        """Compute control action for current state x with optional reference.

        If x_ref is provided, the controller works on the error x_tilde = x - x_ref
        and returns u = -K0 x_tilde + u_ref. If u_ref is None and x_ref is given,
        the method will attempt to compute a steady-state input u_ref that makes
        x_ref an equilibrium (u_ref = pinv(B) (I-A) x_ref). If neither x_ref nor
        u_ref are given, this reduces to the original regulation law u = -K0 x.
        """
        x = np.atleast_1d(x)
        K0 = self.Ks[0]

        if x_ref is None and u_ref is None:
            # Original behavior
            return -K0 @ x

        # ensure proper shapes
        if x_ref is not None:
            x_ref = np.atleast_1d(x_ref)
            x_tilde = x - x_ref
        else:
            x_tilde = x

        if u_ref is None and x_ref is not None:
            # compute steady-state input that makes x_ref an equilibrium if possible
            # Solve (I - A) x_ref = B u_ref  =>  u_ref = pinv(B) (I - A) x_ref
            try:
                I_minus_A = np.eye(self.A.shape[0]) - self.A
                b = I_minus_A @ x_ref
                # Regularized least-squares: solve (B^T B + reg I) u = B^T b
                m = self.B.shape[1]
                BtB = self.B.T @ self.B
                # small regularization scaled to BtB norm for numerical stability
                reg = 1e-8 * max(1.0, np.linalg.norm(BtB))
                try:
                    u_ref = np.linalg.solve(BtB + reg * np.eye(m), self.B.T @ b)
                except Exception:
                    # fallback to pseudoinverse if solve fails
                    u_ref = np.linalg.pinv(self.B) @ b
                # If the computed feedforward is unreasonably large, discard it and
                # fallback to pure feedback control to avoid destabilizing the loop.
                max_u_norm = 1e3
                if np.linalg.norm(u_ref) > max_u_norm:
                    print("Warning: computed steady-state input u_ref is very large; ignoring u_ref and using feedback only.")
                    u_ref = None
            except Exception:
                u_ref = None

        if u_ref is None:
            u_ref = np.zeros(self.B.shape[1])

        u = -K0 @ x_tilde + u_ref
        return u

    def simulate(self, x0: np.ndarray, steps: int, disturbances: Optional[np.ndarray] = None,
                 x_target: Optional[np.ndarray] = None, u_ref: Optional[np.ndarray] = None,
                 verbose: bool = False,
                 stop_on_divergence: bool = True,
                 divergence_factor: float = 10.0,
                 divergence_patience: int = 3,
                 save_path: Optional[str] = None,
                 return_history: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """Simulate the closed-loop system starting from x0 for `steps` steps.

        disturbances: optional array of shape (steps, n) added to the state update.
        x_target: optional desired target state (constant) to track. If provided,
            the controller will compute u_t = -K0 (x - x_target) + u_ref where u_ref
            is either supplied or computed as a steady-state input that makes
            x_target an equilibrium.

        Returns (Xs, Us) where Xs has shape (steps+1, n) and Us has shape (steps, m).
        """
        x = np.atleast_1d(x0).copy()
        n = x.shape[0]
        m = self.B.shape[1]

        Xs = np.zeros((steps + 1, n))
        Us = np.zeros((steps, m))
        errs = np.zeros((steps + 1,))
        Xs[0] = x

        # initialize divergence detection
        if x_target is not None:
            init_err = np.linalg.norm(x - x_target)
        else:
            init_err = np.linalg.norm(x)
        if init_err <= 0:
            init_err = 1e-12
        prev_err = init_err
        increases = 0
        errs[0] = init_err

        for k in range(steps):
            # compute and optionally print error (to target or to zero)
            if x_target is not None:
                err = np.linalg.norm(x - x_target)
            else:
                err = np.linalg.norm(x)
            if verbose:
                print(f"Step {k:3d}: error norm = {err}")

            # divergence checks
            if stop_on_divergence:
                # rapid absolute divergence relative to initial error
                if err > divergence_factor * init_err:
                    if verbose:
                        print(f"Early stop: error {err:.4e} exceeded divergence_factor {divergence_factor} * init_err {init_err:.4e}")
                    Xs_slice = Xs[:k+1].copy()
                    Us_slice = Us[:k].copy()
                    errs_slice = errs[:k+1].copy()
                    if save_path is not None:
                        np.savez_compressed(save_path, Xs=Xs_slice, Us=Us_slice, errs=errs_slice)
                        if verbose:
                            print(f"Saved simulation history to {save_path}")
                    if return_history:
                        return Xs_slice, Us_slice, errs_slice
                    return Xs_slice, Us_slice

                # monotonic growth detection (patience)
                if err > prev_err:
                    increases += 1
                else:
                    increases = 0
                if increases >= divergence_patience and err > 2.0 * init_err:
                    if verbose:
                        print(f"Early stop: error increased for {increases} consecutive steps and reached {err:.4e}")
                    Xs_slice = Xs[:k+1].copy()
                    Us_slice = Us[:k].copy()
                    errs_slice = errs[:k+1].copy()
                    if save_path is not None:
                        np.savez_compressed(save_path, Xs=Xs_slice, Us=Us_slice, errs=errs_slice)
                        if verbose:
                            print(f"Saved simulation history to {save_path}")
                    if return_history:
                        return Xs_slice, Us_slice, errs_slice
                    return Xs_slice, Us_slice

            prev_err = err

            # choose control with optional reference tracking
            if x_target is not None:
                u = self.control_with_ref(x, x_ref=x_target, u_ref=u_ref)
            else:
                u = self.control(x)

            Us[k] = u
            w = disturbances[k] if (disturbances is not None) else 0.0
            x = self.A @ x + self.B @ u + w
            Xs[k + 1] = x
            # record next-step error
            if x_target is not None:
                next_err = np.linalg.norm(x - x_target)
            else:
                next_err = np.linalg.norm(x)
            errs[k + 1] = next_err

        # save final history if requested
        if save_path is not None:
            np.savez_compressed(save_path, Xs=Xs, Us=Us, errs=errs)
            if verbose:
                print(f"Saved simulation history to {save_path}")

        if return_history:
            return Xs, Us, errs
        return Xs, Us


def make_B_from_A(A: np.ndarray, strategy: str = "top_degree", m: int = 1, scale: float = 0.1, base_dir: Optional[str] = None) -> np.ndarray:
    """Create a B matrix when only A is available.

    Strategies:
      - 'single': place single actuator on node 0 (column with 1 at node0)
      - 'identity': use identity (m = n)
      - 'random': random dense B scaled by `scale`
      - 'top_degree': if degree file exists in workspace (degree_data.csv),
           place actuators on top-m degree nodes; otherwise falls back to 'random'
    """
    n = A.shape[0]
    if strategy == "single":
        B = np.zeros((n, m))
        B[0, 0] = 1.0
        return B
    if strategy == "identity":
        # return first m columns of identity
        return np.eye(n, m)
    if strategy == "random":
        rng = np.random.default_rng(1)
        return scale * rng.standard_normal((n, m))
    if strategy == "top_degree":
        # Try to read degree_data.csv from nearby folders
        candidates = []
        if base_dir is None:
            base_dir = os.path.dirname(__file__)
        candidates.extend([
            os.path.join(base_dir, "o2", "degree_data.csv"),
            os.path.join(base_dir, "output", "degree_data.csv"),
            os.path.join(base_dir, "o2", "degree_data.tsv"),
        ])
        degrees = None
        for c in candidates:
            if os.path.exists(c):
                try:
                    degs = []
                    with open(c, newline='') as fh:
                        reader = csv.reader(fh)
                        for row in reader:
                            if not row:
                                continue
                            try:
                                degs.append(float(row[-1]))
                            except Exception:
                                continue
                    if len(degs) >= n:
                        degrees = np.array(degs[:n])
                        print(f"Loaded degrees from {c}")
                        break
                except Exception:
                    degrees = None
        if degrees is None:
            print("degree file not found or invalid, using random B")
            return make_B_from_A(A, strategy="random", m=m, scale=scale)
        idx = np.argsort(-degrees)[:m]
        B = np.zeros((n, m))
        for j, i_node in enumerate(idx):
            B[i_node, j] = 1.0
        return B
    # default fallback
    return make_B_from_A(A, strategy="random", m=m, scale=scale)


def reconstruct_A_from_npz(path: str) -> Optional[np.ndarray]:
    """Try to reconstruct A from an npz file that may contain a CSR matrix.

    Returns a dense ndarray if successful, otherwise None.
    """
    try:
        data = np.load(path)
        keys = set(data.files)
        if set(['data', 'indices', 'indptr', 'shape']).issubset(keys):
            arr_data = data['data']
            arr_indices = data['indices']
            arr_indptr = data['indptr']
            arr_shape = tuple(int(x) for x in data['shape'])
            nrows = arr_shape[0]
            dense = np.zeros(arr_shape, dtype=arr_data.dtype)
            for i in range(nrows):
                start = int(arr_indptr[i])
                end = int(arr_indptr[i + 1])
                if end > start:
                    cols = arr_indices[start:end]
                    vals = arr_data[start:end]
                    dense[i, cols] = vals
            return dense
        # fallback: find a square 2D array
        for k in data.files:
            a = data[k]
            if getattr(a, 'ndim', None) == 2 and a.shape[0] == a.shape[1]:
                return a
    except Exception:
        return None
    return None


if __name__ == "__main__":
    # Example: try to load A/B from workspace. Fall back to a small random system.
    base = os.path.dirname(__file__)
    example_paths = [
        os.path.join(base, "output", "A_final_stable.npz"),
        os.path.join(base, "o2", "A_final_stable.npz"),
    ]

    A = None
    B = None
    for p in example_paths:
        if os.path.exists(p):
            print('Found candidate:', p)
            # first try reconstructing A from CSR-format npz
            A_candidate = reconstruct_A_from_npz(p)
            if A_candidate is not None:
                A = A_candidate
                print(f"Reconstructed A from {p} (shape {A.shape})")
            else:
                try:
                    data = np.load(p)
                    for key in data.files:
                        if key.lower().startswith('a'):
                            A = data[key]
                            print(f"Loaded A from {p} key {key} (shape {A.shape})")
                        if key.lower().startswith('b'):
                            B = data[key]
                except Exception:
                    pass
            if A is not None:
                break

    if A is None:
        print("Could not load A from workspace files, using random example.")
        n = 6
        rng = np.random.default_rng(0)
        M = rng.standard_normal((n, n)) * 0.1
        A = np.eye(n) + M

    if B is None:
        print("B not found — building B from A")
        n = A.shape[0]
        m = min(max(1, n // 10), 1000)
        B = make_B_from_A(A, strategy="top_degree", m=m, scale=0.1, base_dir=base)
        # save generated B for reproducibility
        save_dir = os.path.join(base, 'output') if os.path.exists(os.path.join(base, 'output')) else os.path.join(base, 'o2')
        os.makedirs(save_dir, exist_ok=True)
        B_path = os.path.join(save_dir, 'B_generated.npy')
        np.save(B_path, B)
        print('Saved generated B to', B_path)

    # Build cost matrices
    n = A.shape[0]
    m = B.shape[1]
    Q = np.eye(n)
    R = 0.1 * np.eye(m)
    Qf = 10.0 * np.eye(n)

    print(f"System size: n={n}, m={m}")
    print(np.unique(B))

    horizon = 20
    mpc = MPCController(A, B, Q, R, N=horizon, Qf=Qf)

    # initial state
    x0 = np.ones(n)*0.5

    # desired target state: try to move activity to node 1 (or last node) instead
    x_target = np.zeros(n)

    # simulate closed-loop while tracking x_target
    Xs, Us = mpc.simulate(x0, steps=50, x_target=x_target, verbose=True, save_path="./simulation_history.npz")
    final = Xs[-1]
    print(f"Simulated {Xs.shape[0]-1} steps. Final state norm: {np.linalg.norm(final):.4e}")
    print(f"Distance to target: {np.linalg.norm(final - x_target)}")
