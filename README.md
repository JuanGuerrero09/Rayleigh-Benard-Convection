# Rayleigh-Benard-convection-

This report provides a comprehensive technical breakdown of **Rayleigh-Bénard Convection (RBC)**, covering its physical principles, the mathematical derivation of the non-dimensional governing equations, and the numerical implementation using the **SIMPLE algorithm**.

---

## 1. Physical Description: What is Rayleigh-Bénard Convection?

Rayleigh-Bénard convection occurs when a fluid is trapped between two horizontal plates, where the bottom plate is hotter than the top. This setup creates a temperature gradient that competes with gravity.

1. **Thermal Expansion:** The fluid at the bottom heats up, becomes less dense, and gains buoyancy.
2. **The Instability:** As long as the temperature difference is small, heat is transferred only by **conduction** (the fluid remains still).
3. **Convective Cells:** Once the buoyancy force overcomes the internal friction (viscosity) and thermal dissipation (diffusion), the fluid breaks into organized patterns known as **Bénard Cells**.

---

## 2. Governing Equations (Dimensional)

We assume the fluid is Newtonian and incompressible. We apply the **Boussinesq Approximation**, which treats density as constant everywhere except in the buoyancy term of the momentum equation.

### Mass Conservation (Continuity)

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0$$

### Momentum (Navier-Stokes)

$$\rho_0 \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} \right) = -\nabla p + \mu \nabla^2 \mathbf{u} + \rho(T)\mathbf{g}$$

Using $\rho(T) = \rho_0 [1 - \alpha(T - T_{ref})]$, the vertical component () becomes:

$$\frac{\partial v}{\partial t} + (\mathbf{u} \cdot \nabla)v = -\frac{1}{\rho_0}\frac{\partial p}{\partial y} + \nu \nabla^2 v + g \alpha (T - T_{cold})$$

### Energy (Scalar Transport)

$$\frac{\partial T}{\partial t} + (\mathbf{u} \cdot \nabla)T = \kappa \nabla^2 T$$

---

## 3. Non-Dimensionalization

To generalize the model for any fluid, we define scales:

* **Length:** $h$  (Distance between plates)
* **Time:** $h^2/\kappa$  (Thermal diffusion time)
* **Velocity:** $\kappa/h$
* **Temperature:** $\theta = (T - T_{cold}) / (T_{hot} - T_{cold})$

Substituting these into the dimensional equations yields the **Dimensionless Governing Equations**:

1. **Continuity:** $\nabla \cdot \mathbf{u} = 0$
2. **Momentum:** $\frac{1}{Pr} \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} \right) = -\nabla p + \nabla^2 \mathbf{u} + Ra \cdot \theta \cdot \mathbf{\hat{k}}$
3. **Energy:**  $\frac{\partial \theta}{\partial t} + (\mathbf{u} \cdot \nabla)\theta = \nabla^2 \theta$

### Key Dimensionless Numbers

* **Prandtl Number ($Pr = \nu/\kappa$):** Ratio of momentum diffusivity to thermal diffusivity.
* **Rayleigh Number ($Ra = \frac{g \alpha \Delta T h^3}{\nu \kappa}$):** Determines the strength of the convective drive.

---

## 4. Numerical Implementation: The SIMPLE Algorithm

The **SIMPLE** (Semi-Implicit Method for Pressure-Linked Equations) algorithm is used to solve the coupling between velocity and pressure.

### Step 1: Solving the Scalar (Temperature/Density)

The energy equation is solved first to obtain the temperature field $\theta^{n+1}$. In a discretized 2D grid (indices $i, j$), the advection-diffusion equation is written as a linear system:$$A_P \theta_{i,j} = A_E \theta_{i+1,j} + A_W \theta_{i-1,j} + A_N \theta_{i,j+1} + A_S \theta_{i,j-1} + b_\theta$$

This forms a **Pentadiagonal Matrix** $A_\theta \Theta = B_\theta$, which can be solved using solvers like TDMA or Conjugate Gradient.

### Step 2: The Predictor Step (Momentum)

Using the pressure from the previous time step ($p^n$), we solve the momentum equations to find **intermediate velocities** ($u^*, v^*$):

$$A_P^u u^*_{i,j} = \sum A_{nb} u^*_{nb} + b_u - \left( \frac{p_{i+1,j}^n - p_{i-1,j}^n}{2\Delta x} \right)$$

**Crucial:** These intermediate velocities $u^*, v^*$ do **not** satisfy the continuity equation ($\nabla \cdot \mathbf{u}^* \neq 0$).

### Step 3: Pressure Correction Equation

We define a pressure correction $p'$ such that $p^{n+1} = p^n + p'$. We derive a Poisson Equation for $p'$ by substituting the velocity correction into the continuity equation:

$$\nabla^2 p' = \frac{\nabla \cdot \mathbf{u}^*}{\Delta t}$$

In matrix form: $A_{p'} P' = R_{p}$, where $R_p$ is the "residual" (the mass imbalance of the starred velocities).

### Step 4: Correction and Convergence

1. **Correct Velocity:** $u^{n+1} = u^* - \frac{\Delta t}{\Delta x}(p'_{i+1} - p'_i)$
2. **Correct Pressure:** $p^{n+1} = p^n + \omega p'$ (where $\omega$ is an under-relaxation factor).
3. **Repeat:** If the mass residual is above a threshold, return to Step 2. If it is converged, move to the next time step.

---

## 5. Summary of the Computational Loop

```text
Initialize u, v, p, theta
While t < t_end:
    1. Solve Energy Eq. -> Get theta(n+1)
    2. Solve Momentum Eq. (using p^n) -> Get u*, v*
    3. Solve Pressure Poisson Eq. -> Get p'
    4. Correct u, v and p using p'
    5. Check Convergence (div u < epsilon)
       If No: Repeat from 2
       If Yes: t = t + dt, Update variables

```

This structured approach allows us to simulate how the fluid begins to rotate, forming the iconic convection rolls as  increases beyond the critical value.
