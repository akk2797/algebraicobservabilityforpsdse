# Algebraic Observability for Power System Dynamic State Estimation

MATLAB and Wolfram Mathematica implementations accompanying the paper **“Algebraic Observability for Power System Dynamic State Estimation.”**

The MATLAB scripts evaluate the Jacobian rank of a polynomial differential-algebraic equation (DAE) formulation using successive derivatives of the stator equations. The Mathematica notebooks investigate Smith-form properties of differential-operator matrices for the single-machine DAE and its lifted representation.

## Requirements

- MATLAB with the **Symbolic Math Toolbox**.
- **MATPOWER `case39.m`** on the MATLAB path for the IEEE 39-bus example. The full MATPOWER installation is not required by the supplied IEEE 39-bus script.
- Wolfram Mathematica for the `.nb` notebooks.

## MATLAB files

| File | System and measurements | Main result |
| --- | --- | --- |
| `case1asinglemachine.m` | Single seven-state machine; terminal voltage and current phasors measured | Full symbolic column rank **17/17** at prolongation order **q = 3** |
| `case2a_wscc9_symbolic.m` | WSCC 9-bus, three machines; PMUs at generator buses **1, 2, 3** | Full symbolic column rank **69/69** at **q = 3** |
| `case2b_wscc9_symbolic.m` | WSCC 9-bus, three machines; PMUs at generator buses **1, 2** only | Full column rank **not reached** through **q = 3** in the supplied script |
| `case3_ieee39_symbolic.m` | IEEE 39-bus, ten machines; PMUs at generator buses **30–39** | Full symbolic column rank **257/257** at **q = 3** |

All MATLAB scripts use the original network equations when a network is present; **no Kron reduction** is performed. Higher machine derivatives arising in the stator prolongations are replaced using the machine differential equations.

### Running the MATLAB scripts

Open the desired `.m` file and choose **Run**, or execute its filename in the MATLAB Command Window. For the IEEE 39-bus case, first confirm that `case39.m` is available:

```matlab
which case39 -all
mpc = case39;
[size(mpc.bus,1), size(mpc.gen,1), size(mpc.branch,1)]
```

The expected network counts are `[39 10 46]`. If `which case39 -all` shows several copies, MATLAB uses the first non-shadowed file.

The rank tests use MATLAB's symbolic `jacobian(P, Theta)` and `rank(J)`. No numerical SVD threshold is used for these tests.

### Expected rank progression

| Case | q = 0 | q = 1 | q = 2 | q = 3 |
| --- | --- | --- | --- | --- |
| Single machine (Case 1A) | 12/17 | 14/17 | 16/17 | **17/17** |
| WSCC, three PMUs (Case 2A) | 54/69 | 60/69 | 66/69 | **69/69** |
| WSCC, two PMUs (Case 2B) | 60/75 | 68/81 | 76/87 | 82/93 |
| IEEE 39-bus (Case 3) | 207/257 | 227/257 | 247/257 | **257/257** |

Entries show **Jacobian rank / number of latent indeterminates**. In Case 2B, the unknown count increases with each prolongation because the terminal voltage and current phasors at generator bus 3 are unmeasured. In an additional run extended to **q = 4**, the result was **86/99**; the supplied `case2b_wscc9_symbolic.m` has `qMax = 3`, so change that setting if reproducing the extra order.

A full generic symbolic column rank supports algebraic observability under the assumptions of the selected DAE model and measurements. Failure to reach full rank by a finite tested order does **not** establish non-observability for all higher orders.

## Mathematica notebooks — Case 1B

| File | Formulation | Operator-matrix size |
| --- | --- | --- |
| `Case1B_Original_DAE_clean.nb` | Seven-state machine with original trigonometric stator equations and explicit terminal measurement identities | **13 × 11** |
| `Case1B_Lifted_DAE_clean.nb` | Polynomial rotor-angle lift, with sine/cosine differential equations and unit-circle constraint | **17 × 14** |

Open a `.nb` file in Mathematica and select **Evaluation → Evaluate Notebook**. Each notebook constructs the Jacobians of its equations with respect to variables and their time derivatives, forms **P(s)**, and applies `PolynomialSmithDecomposition` after substituting a specified exact DAE-consistent equilibrium.

The notebooks report the dimensions of **P(s)** and **S(s)**, the Smith form and its diagonal entries, the check **U.P(s).V = S(s)**, and whether all diagonal entries are units. Nonzero constant diagonal entries may be normalized to 1 by invertible polynomial row operations; large integer constants are not observability ranks.

**Scope of the Mathematica result:** The Smith decomposition is for an **equilibrium-specialized polynomial matrix in a commuting variable `s`**. It should not be described as a full noncommutative Jacobson-form or time-varying hyper-regularity proof without additional analysis. The dimensions stated above describe the implemented measurement-augmented matrices and should be distinguished from any differently sized matrix printed in the paper.

## Modeling notes

- Machines are represented by a fourth-order synchronous generator model with an IEEE Type-1 exciter (seven original differential equations per machine).
- Rotor-angle trigonometric functions are polynomially lifted using sine/cosine variables and unit-circle relations.
- Terminal PMU phasors and known inputs are treated as known quantities; unmeasured network variables are included among the latent variables.
- The present scripts treat the excitation saturation coefficient as specified in the current model. An exponential-saturation lifting is **not** included.
- Numerical machine parameters for a physical benchmark should not be inferred from the symbolic rank result alone.

## Reproducibility

Run each script from a clean session and retain the Command Window output with the source files. The MATLAB rank progressions above were reproduced in the submitted implementations; the Mathematica notebooks should be run and their decomposition checks inspected before using their outputs as supporting evidence.
