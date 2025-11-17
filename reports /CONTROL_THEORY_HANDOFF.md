# Control Theory Team - Handoff Document

## Quick Start Guide for Control Design

**System:** Gene Regulatory Network (GRN) for Breast Cancer  
**State Dimension:** n = 8,378 genes  
**System Matrix:** A_final_stable.npz (stable, sparse, signed)  
**Status:** Ready for control design

---

## What You're Getting

### 1. System Matrix (A)

```python
# Load the matrix
from scipy.sparse import load_npz
A = load_npz("o2/A_final_stable.npz")

# Properties:
Shape: 8378 × 8378
Non-zeros: 25,134 (0.036% dense)
Format: CSR sparse matrix
Data type: float32
```

**Key Properties:**
- ✅ **Stable:** All eigenvalues have Re(λ) < 0
- ✅ **Sparse:** 99.96% zeros (efficient computation)
- ✅ **Signed:** Positive (activation) + Negative (repression)
- ✅ **Validated:** Spectral abscissa = -1.16

---

## System Formulation

### Continuous-Time Linear System

```
dx/dt = A·x(t) + B·u(t)

State equation:
  x(t) ∈ ℝ^8378 : Gene expression levels
  A ∈ ℝ^8378×8378 : Regulatory interactions (provided)
  B ∈ ℝ^8378×m : Control input matrix (you design)
  u(t) ∈ ℝ^m : Control inputs (drug dosages)
```

### What's Provided vs What You Design

| Component | Status | Your Task |
|-----------|--------|-----------|
| **A matrix** | ✅ Provided | Use as-is |
| **Gene names** | ✅ Provided | Map indices to biology |
| **B matrix** | ❌ Not designed | **Design this** |
| **Control law u(t)** | ❌ Not designed | **Design this** |
| **Target state x_ref** | ❌ Not defined | **Define this** |

---

## Your Tasks

### Task 1: Define Control Objective

**Options:**

**A. Stabilization (Already Done)**
- System is already stable
- No control needed for stability
- But can improve convergence rate

**B. Set-Point Regulation**
```
Goal: Drive x(t) → x_ref

Example:
  x_ref = healthy gene expression profile
  
Control law:
  u(t) = -K(x(t) - x_ref)
```

**C. Trajectory Tracking**
```
Goal: Follow time-varying trajectory x_d(t)

Example:
  x_d(t) = gradual transition from cancer to healthy
  
Control law:
  u(t) = -K(x(t) - x_d(t)) + u_ff(t)
```

**D. Optimal Control**
```
Goal: Minimize cost function

Example:
  J = ∫[x'Qx + u'Ru] dt
  
  Q: State cost (penalize deviation from equilibrium)
  R: Control cost (penalize large drug dosages)
  
Solution: LQR (Linear Quadratic Regulator)
```

---

### Task 2: Design B Matrix (Control Input Matrix)

**Question:** Which genes can we control (perturb with drugs)?

**Strategy 1: Single-Input Control (m=1)**
```python
# Select one master regulator
control_gene = "TP53"  # Example
gene_idx = genes.index(control_gene)

B = np.zeros((8378, 1))
B[gene_idx, 0] = 1.0

# Interpretation: We can directly control TP53 expression
```

**Strategy 2: Multi-Input Control (m=10)**
```python
# Select top-10 hub genes (master regulators)
hub_genes = ["TP53", "MYC", "BRCA1", "ESR1", ...]
m = len(hub_genes)

B = np.zeros((8378, m))
for i, gene in enumerate(hub_genes):
    gene_idx = genes.index(gene)
    B[gene_idx, i] = 1.0

# Interpretation: We can control 10 key genes simultaneously
```

**Strategy 3: Minimum Dominating Set**
```python
# Find smallest set of genes that can reach all others
# Use graph algorithms (see graph_stats.json)
# Driver nodes: 7,899 candidates (94% of network)
```

**How to Choose:**
- **Biological:** Which genes have drug targets?
- **Structural:** Which genes have high out-degree (hubs)?
- **Controllability:** Which genes maximize rank(C)?

---

### Task 3: Check Controllability

**Controllability Matrix:**
```python
from scipy.linalg import matrix_rank

# Build controllability matrix
C = [B]
AB = A @ B
for i in range(1, n):  # n = 8378 (expensive!)
    AB = A @ AB
    C.append(AB)

C = np.hstack(C)

# Check rank
rank_C = matrix_rank(C)
is_controllable = (rank_C == n)

print(f"Rank(C) = {rank_C}")
print(f"Controllable: {is_controllable}")
```

**Warning:** This is computationally expensive for n=8,378!

**Alternative: Structural Controllability**
```python
# Use graph theory (faster)
# Check if selected genes form a dominating set
# See: Liu et al. (2011) Nature paper
```

---

### Task 4: Design Control Law

**Option A: State Feedback (LQR)**

```python
from scipy.linalg import solve_continuous_are

# Define cost matrices
Q = np.eye(n)  # State cost (identity)
R = np.eye(m)  # Control cost (identity)

# Solve Riccati equation: A'P + PA - PBR^(-1)B'P + Q = 0
P = solve_continuous_are(A.toarray(), B, Q, R)

# Compute feedback gain: K = R^(-1)B'P
K = np.linalg.inv(R) @ B.T @ P

# Control law
def control_law(x):
    return -K @ x

# Closed-loop system: dx/dt = (A - BK)x
A_cl = A - B @ K
```

**Option B: Model Predictive Control (MPC)**

```python
# Solve optimization at each time step:
# min  Σ[x'Qx + u'Ru]
# s.t. x(k+1) = Ax(k) + Bu(k)
#      u_min ≤ u(k) ≤ u_max
#      x_min ≤ x(k) ≤ x_max

# Use CVXPY or similar optimization library
```

---

### Task 5: Simulate Controlled System

```python
from scipy.integrate import odeint

# Define closed-loop dynamics
def closed_loop_dynamics(x, t, A, B, K):
    u = -K @ x  # Control law
    dxdt = A @ x + B @ u
    return dxdt

# Initial condition (perturbed from equilibrium)
x0 = np.random.randn(n) * 0.1

# Time span
t = np.linspace(0, 10, 1000)

# Simulate
x_traj = odeint(closed_loop_dynamics, x0, t, args=(A, B, K))

# Plot
import matplotlib.pyplot as plt
plt.plot(t, np.linalg.norm(x_traj, axis=1))
plt.xlabel("Time")
plt.ylabel("||x(t)||")
plt.title("Controlled System Response")
plt.show()
```

---

## Example: Cancer Therapy Control

### Objective
Suppress oncogene (MYC) while maintaining tumor suppressor (TP53)

### Target State
```python
x_ref = np.zeros(n)
x_ref[genes.index("TP53")] = 1.0   # High expression
x_ref[genes.index("MYC")] = -1.0   # Low expression
```

### Control Genes
```python
# Select genes that regulate TP53 and MYC
control_genes = ["ESR1", "BRCA1", "ATM"]  # Example
```

### Control Law
```python
# LQR with set-point tracking
def control_law(x, x_ref, K):
    return -K @ (x - x_ref)
```

### Expected Outcome
```
t=0:   MYC high, TP53 low (cancer state)
t=5:   MYC decreasing, TP53 increasing (transition)
t=10:  MYC low, TP53 high (healthy state)
```

---

## Practical Considerations

### 1. Dimensionality Reduction

**Problem:** 8,378 states is too large for some control methods

**Solution:** Model reduction
```python
from sklearn.decomposition import PCA

# Reduce to k principal components
pca = PCA(n_components=100)
x_reduced = pca.fit_transform(x_traj)

# Design control in reduced space
# Then map back to full space
```

### 2. Sparse Actuation

**Problem:** Can't control all 8,378 genes (limited drug targets)

**Solution:** Select m << n control genes
```python
# Use biological knowledge
# - Known drug targets
# - Master regulators (high out-degree)
# - Driver nodes (no incoming edges)
```

### 3. Constraints

**Problem:** Gene expression has physical bounds

**Solution:** Add constraints to optimization
```python
# State constraints
x_min = -3.0  # Min expression (log-scale)
x_max = +3.0  # Max expression (log-scale)

# Control constraints
u_min = 0.0   # Min drug dosage
u_max = 1.0   # Max drug dosage
```

### 4. Uncertainty

**Problem:** Model mismatch, noise, disturbances

**Solution:** Robust control
```python
# H-infinity control
# Robust MPC
# Adaptive control
```

---

## Deliverables

### What We Expect

1. **B Matrix Design**
   - Which genes selected as control inputs?
   - Justification (biological + structural)

2. **Controllability Analysis**
   - Is system controllable with chosen B?
   - Rank of controllability matrix

3. **Control Law**
   - Feedback gain K (for LQR)
   - Or MPC formulation

4. **Simulations**
   - Closed-loop response
   - Convergence to target state
   - Control effort (u(t) trajectory)

5. **Analysis**
   - Stability of closed-loop system
   - Performance metrics (settling time, overshoot)
   - Robustness analysis

---

## Code Templates

### Template 1: Load System

```python
import numpy as np
import pandas as pd
from scipy.sparse import load_npz

# Load system matrix
A = load_npz("o2/A_final_stable.npz")
genes = pd.read_csv("o2/genes_final.csv")["gene"].tolist()

print(f"System size: {A.shape[0]} genes")
print(f"System is stable: {True}")  # Pre-verified
```

### Template 2: Design B Matrix

```python
# Select control genes
control_genes = ["TP53", "MYC", "BRCA1"]  # Example
m = len(control_genes)

# Build B matrix
B = np.zeros((A.shape[0], m))
for i, gene in enumerate(control_genes):
    if gene in genes:
        gene_idx = genes.index(gene)
        B[gene_idx, i] = 1.0
    else:
        print(f"Warning: {gene} not in gene list")

print(f"B matrix shape: {B.shape}")
print(f"Control inputs: {m}")
```

### Template 3: LQR Design

```python
from scipy.linalg import solve_continuous_are

# Cost matrices
Q = np.eye(A.shape[0])
R = np.eye(m)

# Solve Riccati equation
P = solve_continuous_are(A.toarray(), B, Q, R)

# Feedback gain
K = np.linalg.inv(R) @ B.T @ P

print(f"Feedback gain shape: {K.shape}")
```

### Template 4: Simulate

```python
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Closed-loop dynamics
def dynamics(x, t):
    u = -K @ x
    dxdt = A @ x + B @ u
    return dxdt

# Initial condition
x0 = np.random.randn(A.shape[0]) * 0.1

# Simulate
t = np.linspace(0, 10, 1000)
x_traj = odeint(dynamics, x0, t)

# Plot
plt.plot(t, np.linalg.norm(x_traj, axis=1))
plt.xlabel("Time")
plt.ylabel("State Norm")
plt.title("Controlled System Response")
plt.grid(True)
plt.show()
```

---

## Helpful Resources

### Files to Use

1. **A_final_stable.npz** - System matrix (essential)
2. **genes_final.csv** - Gene names (essential)
3. **degree_data.csv** - Hub genes (for B matrix design)
4. **graph_stats.json** - Network properties (for controllability)

### Key Statistics

```json
{
  "n_nodes": 8378,
  "n_edges": 25134,
  "spectral_abscissa": -1.16,
  "n_driver_nodes": 7899,
  "driver_nodes_fraction": 0.94
}
```

### Suggested Control Genes (Hubs)

Based on biological knowledge + network structure:
- TP53 (tumor suppressor, master regulator)
- MYC (oncogene, high out-degree)
- BRCA1 (tumor suppressor, DNA repair)
- ESR1 (estrogen receptor, breast cancer)
- EGFR (growth factor receptor)
- PIK3CA (signaling pathway)
- AKT1 (cell survival)
- PTEN (tumor suppressor)
- RB1 (cell cycle regulator)
- CDKN2A (cell cycle inhibitor)

---

## Questions?

**For system matrix questions:**
- Check README.md (comprehensive documentation)
- Check integration_params.json (pipeline parameters)
- Check stabilization_info.json (eigenvalue details)

**For control design questions:**
- Standard textbooks: Ogata, Franklin et al., Åström & Murray
- LQR: Boyd & Barratt "Linear Controller Design"
- Network control: Liu et al. (2011) Nature paper

**For biological context:**
- TCGA breast cancer data
- Gene function: NCBI Gene database
- Pathways: KEGG, Reactome

---

## Timeline Suggestion

**Week 1:**
- Load system, understand structure
- Select control genes (B matrix)
- Check controllability

**Week 2:**
- Design control law (LQR or MPC)
- Implement simulation
- Validate closed-loop stability

**Week 3:**
- Run simulations
- Analyze performance
- Generate plots

**Week 4:**
- Write report
- Prepare presentation
- Document code

---

## Success Criteria

Your control design is successful if:

1. ✅ **Controllable:** rank(C) = n or structural controllability proven
2. ✅ **Stable:** Closed-loop eigenvalues have Re(λ) < 0
3. ✅ **Converges:** x(t) → x_ref as t → ∞
4. ✅ **Feasible:** Control inputs u(t) are bounded
5. ✅ **Robust:** Works with model uncertainty
6. ✅ **Biologically meaningful:** Control genes are druggable

---

## Contact

For questions about the system matrix or GRN:
- [Your name/email]

For control theory questions:
- [Control theory team lead]

---

**Good luck with the control design!** 🎮🚀

