# G1: Guaranteed lower bound for the exclusion criterion (Theorem G1)

Companion to `DESIGN.md`. That file proves: **if for every s ∈ (0,1) the form
𝒜ₛ is positive on the hyperplane {ℒₛ = 0} ⊂ H¹(K), then λ₁(PSL(2,ℤ[i])\ℍ³) ≥ 1.**
This file turns that infinite-dimensional positivity into a finite, verified
computation. Every constant is explicit and interval-computable; `cr_prototype.py`
mirrors this file lemma-for-lemma (float now, Arb in M3).

Provenance of external inputs:
- (I1) for the Crouzeix–Raviart (CR) interpolation on tetrahedra: proved
  self-contained below (Lemma I1, arbitrary tets, κ_sc = √(1/π²+1/15));
  classical inputs only (Payne–Weinberger 1960, Bebendorf 2003). The sharper
  κ₁ = √(1/π²+1/120) of Carstensen–Puttkammer, *Adaptive guaranteed lower
  eigenvalue bounds with optimal convergence rates*, §2.4 (arXiv:2203.01028;
  sourced there to Carstensen–Zhai–Zhang, SIAM J. Numer. Anal. 58 (2020)
  109–124 and Carstensen–Puttkammer, J. Comput. Math. 38 (2020) 142–175;
  local copies in `refs/`) is available as an option but not load-bearing.
- The GLB template (project through I_CR, Young inequalities, min-structure):
  Carstensen–Gedicke, Math. Comp. 83 (2014) 2605–2629 (`refs/`), and the
  lower-*energy*-bound viewpoint of Carstensen–Dond–Maity–Nataraj, CAMWA 214
  (2026) 27–50, doi:10.1016/j.camwa.2026.03.029. **[TODO]** the CAMWA paper
  is paywalled; PDF to be added to `refs/` by the user; reconcile constants.
- Everything else (trace identity, sliver estimates, pencil reduction) is
  proved here from scratch.

## 0. The sufficient pencil criterion

Notation from DESIGN.md: K = {(x₁,x₂,y) : x₁∈[−½,½], x₂∈[0,½],
y_f(z) ≤ y ≤ Y}, y_f = √(1−|z|²), Y > 1 fixed (Y = 1.25). For v ∈ H¹(K):

    Q(v) = ∫_K |∇v|² y⁻¹ dx dy        (weighted Dirichlet form)
    M(v) = ∫_K v² y⁻³ dx dy           (weighted mass)
    t(v) = ∫_{y=Y} v dx₁dx₂           (Euclidean top-face integral)
    a(v) = ∫_K v y⁻³ dx dy
    ℒ_s(v) = a(v) + κ_c(s) t(v),   κ_c(s) = 1/((1+s)Y²)
    λ = 1−s²,   β(s) = 2(1−s)/Y².

**Lemma P (pencil criterion).** Fix ρ > 0, ε₀ > 0. If for every s ∈ (0,1)

    Q(v) + ρ ℒ_s(v)²  ≥  (1+ε₀) [ λ M(v) + β(s) t(v)² ]      for all v ∈ H¹(K),   (P_s)

then the criterion of DESIGN.md §2 holds (with strict positivity), hence
λ₁ ≥ 1.

*Proof.* Let v ≠ 0 with ℒ_s(v) = 0. Then (P_s) gives
𝒜_s(v) = Q − λM − βt² ≥ ε₀(λM(v) + β t(v)²) ≥ ε₀ λ M(v) > 0. ∎

(P_s) is *unconstrained* and both sides are nonnegative forms — the right
container for nonconforming lower bounds. ρ is a free parameter; correctness
never depends on its value.

## 1. Mesh and geometry

Start from the M0 mapped hexahedral grid (N₁×N₂×N₃ cells over
(ξ₁,ξ₂,ξ₃) ∈ [−½,½]×[0,½]×[0,1], vertical coordinate y = y_f(1−ξ₃)+Yξ₃).
**Floor lift:** every node with ξ₃ = 0 is placed at y = y_f(z) + ℓ instead of
y_f(z), with the lift

    ℓ := ½ · max_cell |D²y_f| · (cell diameter in z)²,  |D²y_f| ≤ y_f⁻³ ≤ 2^{3/2},

evaluated per floor cell (code uses the cell's max of y_f⁻³ and its z-diameter).
Each (possibly non-planar) hexahedral cell is split into 6 tetrahedra by the
Kuhn triangulation with globally consistent diagonals; all tet faces are
planar (each tet is a convex hull of 4 nodes) and the mesh is face-to-face.
K_h := union of the tets; Γ_h := its floor (a piecewise-linear graph
ŷ(z) over the rectangle R = [−½,½]×[0,½]); Σ := K \ K_h = {(z,y) :
y_f(z) ≤ y < ŷ(z)} the sliver; δ(z) := ŷ − y_f ≥ 0, δ̄ := max δ.

**Lemma G (inclusion).** y_f is concave on R. On a floor triangle with
longest z-edge d, the chord interpolant lies below y_f by at most the sag
(max|D²y_f|)·d²/8. With the per-cell lift ℓ := (max_cell y_f⁻³)·d²/8
(|D²y_f| ≤ y_f⁻³) assigned to nodes as the max over adjacent cells,
ŷ ≥ y_f pointwise, i.e. K_h ⊂ K, and δ(z) ≤ max vertex lift =: δ̂ (per face).

*Proof.* Taylor with concavity: chord(z) = Σθᵢy_f(zᵢ) ≥ y_f(z) −
½max|D²y_f| Σθᵢ|zᵢ−z|². The barycentric variance Σθᵢ|zᵢ−z|² at z = Σθᵢzᵢ
equals Σθᵢ|zᵢ|² − |Σθᵢzᵢ|² ≤ d²/4 (maximized at an edge midpoint of the
longest edge). Adding the lifts raises the chord by ≥ min vertex lift ≥ the
sag. Since chord ≤ y_f, also δ = (chord + interp ℓ) − y_f ≤ max vertex
lift. ∎ (Code self-check: ŷ ≥ y_f sampled on every floor face.)

Per-tet data (all computed by code): h_T = diam T, |T|, face areas |F|;
y-range [y_T^-, y_T^+] over T (min/max over the 4 vertices — valid since y is
affine on T); sandwiched weights

    w_T^Q := (y_T^+)⁻¹  (stiffness, from below),
    w_T^M := (y_T^-)⁻³  (mass, from above),

so that Q(v) ≥ Q_pc(v) := Σ_T w_T^Q ‖∇v‖²_{L²(T)} and
M|_{K_h}(v) ≤ M_pc(v) := Σ_T w_T^M ‖v‖²_{L²(T)} for every v ∈ H¹(K_h).

## 2. Tools

**(I1)–(I2) CR interpolation.** V_h := CR¹(𝒯) (piecewise P1, continuous at
face barycenters/means, *no* boundary conditions). The face-mean
interpolation I = I_CR : H¹(K_h) → V_h satisfies, on each tet T, with
e := v − Iv:

    (I2)  ∇(Iv)|_T = Π₀(∇v|_T)   ⇒   ∫_T ∇e·∇w_h = 0 for w_h ∈ V_h,
          hence ‖∇v‖²_T = ‖∇Iv‖²_T + ‖∇e‖²_T  and  Q_pc(v) = Q_pc(Iv) + Q_pc(e);
    (I1)  ‖e‖_{L²(T)} ≤ κ h_T ‖∇e‖_{L²(T)}.

(I2) is immediate from the face-mean property (∫_T ∇e = Σ_F ν_F ∫_F e = 0
and ∇Iv is constant). For (I1) two constants are available:

**Lemma I1 (self-contained, arbitrary tetrahedra).** (I1) holds with

    κ_sc := √(1/π² + 1/15) ≤ 0.40988.

*Proof.* e = v − I_CR v has ∫_F e dσ = 0 on every face F of T. Split
e = (e − ē_T) + ē_T, orthogonal in L²(T): ‖e‖² = ‖e−ē_T‖² + |T| ē_T².
Payne–Weinberger for the convex domain T [Payne–Weinberger, Arch. Rational
Mech. Anal. 5 (1960) 286–292; corrected proof: Bebendorf, Z. Anal. Anwend.
22 (2003) 751–756] gives ‖e−ē_T‖ ≤ (h_T/π)‖∇e‖. For the mean, pick a
vertex P with opposite face F: div((x−P)e) = 3e + (x−P)·∇e, and
∫_{∂T}(x−P)·ν e dσ = (3|T|/|F|)∫_F e dσ = 0 as in Lemma T; hence
3|T| ē_T = −∫_T (x−P)·∇e; Cauchy–Schwarz gives
|T|ē_T² = (1/(9|T|))(∫_T(x−P)·∇e)² ≤ (1/9)(|T|⁻¹∫_T|x−P|²)‖∇e‖² = (m₂/9)‖∇e‖² with
m₂ := |T|⁻¹∫_T|x−P|². Writing x−P = Σ_{i≠P} λᵢ(vᵢ−P) and using
∫_T λᵢλⱼ = |T|(1+δᵢⱼ)/20: m₂ ≤ h_T²(6·1 + 3·2)/20 = (3/5)h_T². So
|T|ē_T² ≤ h_T²/15 · ‖∇e‖², and κ² ≤ 1/π² + 1/15. ∎

**Sharper published constant.** κ₁ := √(1/π² + 1/120) ≤ 0.33115 per
Carstensen–Puttkammer [CP22 §2.4], sourced to Carstensen–Zhai–Zhang, SIAM
J. Numer. Anal. 58 (2020) 109–124. CZZ20 is paywalled and its hypotheses
were not independently audited here, so the *certificate runs with κ_sc by
default* and uses κ₁ only as an optional sharpening. Nothing below depends
on which is chosen; κ enters only scalar constants (γ, τ_h, α_h, σ_h),
never the matrices.

**Lemma T (trace identity).** Let T be a tet, F a face, P the opposite
vertex. For any w ∈ H¹(T):

    ∫_F w² dσ = (|F|/|T|) ∫_T [ w² + (2/3) w (x−P)·∇w ] dx
              ≤ (|F|/|T|) ( ‖w‖²_T + (2/3) h_T ‖w‖_T ‖∇w‖_T ).

*Proof.* div((x−P)w²) = 3w² + 2w(x−P)·∇w; integrate over T. On the three
faces containing P, (x−P)·ν = 0; on F, (x−P)·ν = dist(P, aff F) = 3|T|/|F|.
Cauchy–Schwarz and |x−P| ≤ h_T. ∎

**Lemma E (functional errors).** For any v ∈ H¹(K_h), e = v − I_CR v:

    (E-M)  M_pc(e) ≤ γ² Q_pc(e),               γ := κ · max_T h_T √(w_T^M / w_T^Q);
    (E-t)  |t(e)| ≤ τ_h √(Q_pc(e)),  for any slab height 0 < H_t ≤ Y − max ŷ:
           τ_h := √(|R|/H_t) · κ · max_{T∩slab} (h_T/√w_T^Q)
                  + √(H_t |R| Y / 3),           |R| = ½, slab = {y ≥ Y − H_t};
    (E-a)  |a(e)| ≤ α_h √(Q_pc(e)),            α_h := κ ( Σ_T (w_T^M)² |T| h_T² / w_T^Q )^{1/2};
    (E-ℒ)  |ℒ_{s₀}(e)| ≤ σ_h √(Q_pc(e)),       σ_h := α_h + κ_c(s₀) τ_h.

*Proof.* (E-M): per tet, w_T^M‖e‖² ≤ w_T^M κ²h_T²‖∇e‖² =
(w_T^M κ²h_T²/w_T^Q)·w_T^Q‖∇e‖²; take the max ratio.

(E-t), slab-mean identity: for a.e. z ∈ R and w ∈ H¹,
w(z,Y) = (1/H_t)∫_{Y−H_t}^{Y} w dy + (1/H_t)∫_{Y−H_t}^{Y} (y−Y+H_t) ∂_y w dy
(integrate the second term by parts to check). The slab {Y−H_t ≤ y ≤ Y} ⊂ K_h
provided H_t ≤ Y − max ŷ. Hence with w = e,

    |t(e)| ≤ (1/H_t)‖1‖_{L²(slab)}‖e‖_{L²(slab)}
             + (1/H_t)‖y−Y+H_t‖_{L²(slab)}‖∂_y e‖_{L²(slab)}
           ≤ √(|R|/H_t)·‖e‖_{euc,slab} + √(H_t|R|/3)·‖∇e‖_{euc,slab},

and per tet ‖e‖²_{euc,T} ≤ κ²h_T²‖∇e‖²_{euc,T} = (κ²h_T²/w_T^Q)·w_T^Q‖∇e‖²,
‖∇e‖²_{euc,T} ≤ (1/w_T^Q)w_T^Q‖∇e‖² ≤ Y·w_T^Q‖∇e‖² (w_T^Q ≥ 1/Y). Sum over
the tets meeting the slab and take maxima. (E-a): per tet
|∫_T e y⁻³| ≤ w_T^M |T|^{1/2}‖e‖_T ≤ w_T^M |T|^{1/2}κh_T‖∇e‖_T; Cauchy–Schwarz.
(E-ℒ): triangle inequality. ∎

(The earlier per-face route via Lemma T gives τ_h = O(√h) as well but with a
~2× worse constant; the slab form is what the code uses, optimizing H_t
numerically — any H_t is admissible, so the optimization is correctness-free.)

**Lemma S (sliver via first-layer column means).** Fix a column height
H_s ∈ (0, Y − max ŷ]. For a floor face F let π(F) ⊂ R be its z-projection,
δ̂_F the max vertex lift among its vertices (≥ δ on π(F), Lemma G),
y⁻_F := min_{π(F)} y_f, y⁺_F := max_{π(F)} ŷ + H_s. The columns
col(F) := {(z,y) : z ∈ π(F), ŷ(z) ≤ y ≤ ŷ(z) + H_s} ⊂ K_h are disjoint.
Then for every v ∈ H¹(K), with Q_Σ(v) := ∫_Σ |∇v|² y⁻¹:

    (S3)  M_Σ(v) := ∫_Σ v² y⁻³ ≤ S_M M_{K_h}(v) + S_Q Q_{K_h}(v) + S_Σ Q_Σ(v),
          S_M := max_F 2 (δ̂_F/H_s) (y⁺_F / y⁻_F)³,
          S_Q := max_F 2 δ̂_F (δ̂_F + H_s) y⁺_F / (y⁻_F)³,
          S_Σ := max_F 2 δ̂_F (δ̂_F + H_s) (max_{π(F)} ŷ) / (y⁻_F)³;
    (S4)  |a_Σ(v)| := |∫_Σ v y⁻³| ≤ √(V_Σ) √(M_Σ(v)),   V_Σ := ∫_Σ y⁻³ ≤ ½ δ̄ y_m⁻³,
          y_m := 1/√2, δ̄ := max_F δ̂_F.

*Proof.* Fix z ∈ π(F) and y ∈ [y_f(z), ŷ(z)]. The mean-value identity over
the column fiber [ŷ, ŷ+H_s] and Cauchy–Schwarz give

    v(z,y)² ≤ 2(1/H_s)∫_{col fiber} v² dy′ + 2(δ̂_F+H_s)∫_{[y, ŷ+H_s]} (∂_y v)² dy′,

since |v(z,y)| ≤ (1/H_s)∫_{col}|v| + ∫_y^{ŷ+H_s}|∂_y v| and the second
interval has length ≤ δ̂_F + H_s. Integrating y over the sliver fiber
(length ≤ δ̂_F), then z over π(F):

    ∫_{Σ∩col⁻} v²_{euc} ≤ 2(δ̂_F/H_s)∫_{col(F)} v²_{euc}
                          + 2δ̂_F(δ̂_F+H_s)[∫_{col(F)} + ∫_{Σ_F}] (∂_y v)²_{euc}.

Weight conversions on each region: on the sliver, y⁻³ ≤ (y⁻_F)⁻³; on col(F),
v²_{euc} ≤ (y⁺_F)³ v² y⁻³ and (∂_y v)²_{euc} ≤ y⁺_F |∇v|² y⁻¹; on Σ_F,
(∂_y v)²_{euc} ≤ (max ŷ)|∇v|²y⁻¹. Sum over F (columns disjoint, sliver
pieces disjoint) and take the max coefficients. (S4): Cauchy–Schwarz,
|Σ| ≤ ½δ̄. ∎

## 3. Master theorem

Fix a window W = [s⁻, s⁺] ⊆ [0,1] and a reference s₀ ∈ W. Worst-case data:

    λ⁺ := 1 − (s⁻)²,   β⁺ := 2(1−s⁻)/Y²,   Δκ := max(|κ_c(s⁻)−κ_c(s₀)|, |κ_c(s⁺)−κ_c(s₀)|).

Choose Young parameters θ, θ₂, θ', α ∈ (0,1) and c > 0, a target ν* = 1+ε₀,
and ρ > 0. Define, in this order (all computable from mesh data):

    ω    := 2ρ(1/θ₂ − 1) V_Σ
    (i)   sliver absorption:
          c_Q := 1 − (ω + ν*λ⁺) S_Q          [numerator Q_{K_h} coefficient]
          c_Σ := 1 − (ω + ν*λ⁺) S_Σ          [must be ≥ 0; then Q_Σ dropped]
          λ̃  := ν*λ⁺(1 + S_M) + ω S_M       [effective mass coefficient]
          β̃  := ν*β⁺ + 2ρ(1/θ₂ − 1)Δκ²      [effective trace coefficient]
          ρ̃  := ρ(1 − θ₂)
    (ii)  CR reduction:
          c_e := c_Q − ρ̃(1/θ − 1)σ_h²
          d_e := λ̃(1 + 1/α)γ² + β̃(1 + 1/θ')τ_h²

Discrete check: with the CR matrices Q_pc, M_pc, vectors ℓ₀ (for
ℒ_{s₀}(v_h) = a(v_h) + κ_c(s₀)t(v_h), a with exact weight y⁻³) and t,

    N_h := c_Q Q_pc + ρ̃(1−θ) ℓ₀ℓ₀ᵀ,
    D_h := λ̃(1+α) M_pc + β̃(1+θ') ttᵀ.

**Theorem G1.** Suppose c_Σ ≥ 0, c_e > 0, and

    (a)  N_h − D_h ⪰ 0   (positive semidefinite on the CR space),
    (b)  c_e / d_e ≥ 1.

Then (P_s) holds with 1+ε₀ = ν* for every s ∈ W. If (a),(b) hold for a
family of windows covering (0,1), then by Lemma P, λ₁(PSL(2,ℤ[i])\ℍ³) ≥ 1.

*Proof.* Let v ∈ H¹(K), s ∈ W. Write N := Q(v) + ρℒ_s(v)²,
D := λM(v) + β(s)t(v)². Since λ ≤ λ⁺ and β(s) ≤ β⁺, D ≤ λ⁺M(v) + β⁺t(v)².

Step 1 (window shift + sliver in ℒ). ℒ_s(v) = ℒ_{s₀}(v)|_{K_h} + a_Σ(v) +
(κ_c(s) − κ_c(s₀)) t(v), so by (x+y)² ≥ (1−θ₂)x² − (1/θ₂−1)y² and
(u+w)² ≤ 2u² + 2w², Lemma S (S4):

    ρℒ_s(v)² ≥ ρ̃ ℒ_{s₀,K_h}(v)² − 2ρ(1/θ₂−1)[ V_Σ M_Σ(v) + Δκ² t(v)² ].

Step 2 (sliver in M). M(v) = M_{K_h}(v) + M_Σ(v) and (S3). Collecting, with
Q(v) = Q_{K_h}(v) + Q_Σ(v):

    N − ν*D ≥ c_Q Q_{K_h}(v) + c_Σ Q_Σ(v) + ρ̃ ℒ_{s₀,K_h}(v)²
              − λ̃ M_{K_h}(v) − β̃ t(v)².

c_Σ ≥ 0 lets us drop Q_Σ. It remains to show, for all v ∈ H¹(K_h),

    c_Q Q_{K_h}(v) + ρ̃ ℒ₀(v)² ≥ λ̃ M_{K_h}(v) + β̃ t(v)².       (★)

Step 3 (pc sandwich). Q_{K_h} ≥ Q_pc and M_{K_h} ≤ M_pc, so it suffices to
prove (★) with Q_pc, M_pc.

Step 4 (CR projection). Split v = Iv + e. By (I2), Q_pc(v) = Q_pc(Iv) + Q_pc(e).
By Lemma E and Young:

    ℒ₀(v)² ≥ (1−θ)ℒ₀(Iv)² − (1/θ−1)σ_h² Q_pc(e),
    M_pc(v) ≤ (1+α)M_pc(Iv) + (1+1/α)γ² Q_pc(e),
    t(v)²   ≤ (1+θ')t(Iv)² + (1+1/θ')τ_h² Q_pc(e).

Therefore

    LHS(★) ≥ [c_Q Q_pc(Iv) + ρ̃(1−θ)ℒ₀(Iv)²] + c_e Q_pc(e),
    RHS(★) ≤ [λ̃(1+α)M_pc(Iv) + β̃(1+θ')t(Iv)²] + d_e Q_pc(e).

By (a), the bracketed terms satisfy LHS-bracket ≥ RHS-bracket (Iv ∈ V_h);
by (b), c_e Q_pc(e) ≥ d_e Q_pc(e). Adding proves (★). ∎

Remarks.
- **What the code must certify** per window: the scalar inequalities
  c_Σ ≥ 0, c_e/d_e ≥ 1, and the matrix inequality N_h − D_h ⪰ 0 (in M3 via
  interval Cholesky of N_h − D_h + small shift; in the float prototype via
  Cholesky / smallest-eigenvalue diagnostics). No eigenvalue *solve* is part
  of the certificate — only positive-definiteness checks.
- All quadratic-form matrix entries are per-tet polynomial integrals with
  constant (sandwiched) weights — exactly representable; the only curved-
  weight objects are the *vector* ℓ₀ (entries ∫_T φ_F y⁻³, one rational 3D
  integral per tet, enclosable in Arb) and the scalar geometric constants.
- Asymptotics (sanity): γ, α_h = O(h), τ_h = O(h^{1/2}) (top-face sum has
  O(h⁻²) terms of size O(h³)·h²/(h³) …), δ̄, S_*, V_Σ = O(h²). So
  c_e/d_e = O(h⁻²) → ∞ and N_h−D_h tends to the unperturbed discrete pencil:
  the theorem is not just true but *efficient* — the h-dependence enters only
  through small multiplicative/additive corrections.
- The strictness needed by DESIGN.md §2 comes from ε₀ > 0 in ν* = 1+ε₀
  (Lemma P gives 𝒜_s ≥ ε₀λM > 0 on the hyperplane).

## 4. Constants table (evaluated by `cr_prototype.py` per mesh/window)

| symbol | meaning | formula/source |
|---|---|---|
| κ | CR interpolation constant, 3D tets | default κ_sc = √(1/π²+1/15) ≤ 0.40988 (Lemma I1, self-contained); optional κ₁ = √(1/π²+1/120) [CP22 §2.4, CZZ20] |
| w_T^Q, w_T^M | sandwiched weights | (y_T^+)⁻¹, (y_T^-)⁻³ |
| γ, τ_h, α_h, σ_h | interpolation-error functionals | Lemma E (τ_h: slab form) |
| ℓ, δ̂_F, δ̄ | floor lift, sliver thickness | Lemma G |
| S_M, S_Q, S_Σ, V_Σ | sliver constants | Lemma S |
| λ⁺, β⁺, Δκ | window worst-case | §3 |
| θ, θ₂, θ', α, H_t, H_s, ρ, ν* | free parameters | tuned numerically, correctness-free |
| c_Q, c_Σ, λ̃, β̃, ρ̃, c_e, d_e | derived coefficients | §3 (i)–(ii) |

## 5. M3 verification appendix (implemented in `m3_certify.py`)

**Lemma R (factored interval rank-one).** Let ℓ ∈ ℝⁿ lie componentwise in
[ℓ̂ ± r]. Then for all x: (ℓᵀx)² ≥ (ℓ̂ᵀx)² − 2‖ℓ̂‖₂‖r‖₂‖x‖₂², i.e.
ℓℓᵀ ⪰ ℓ̂ℓ̂ᵀ − 2‖ℓ̂‖‖r‖·I; and (ℓᵀx)² ≤ (ℓ̂ᵀx)² + (2‖ℓ̂‖‖r‖+‖r‖²)‖x‖², i.e.
ℓℓᵀ ⪯ ℓ̂ℓ̂ᵀ + (2‖ℓ̂‖‖r‖+‖r‖²)·I. *Proof:* write ℓ = ℓ̂+δ, |δᵀx| ≤ ‖r‖‖x‖,
|ℓ̂ᵀx| ≤ ‖ℓ̂‖‖x‖. ∎

**ℓ₀ enclosure.** ∫_T φ_a y⁻³ dx by degree-p Taylor of y⁻³ about the tet
mean ȳ: the moments ∫_T φ_a (y−ȳ)^k are exact barycentric-monomial
integrals (∫_T Πλᵢ^{αᵢ} = 6|T|Παᵢ!/(|α|+3)!), and the remainder is the ball
±2|T|·((p+3)!/(2(p+1)!))·y_T^{−(4+p)}·Δ^{p+1}, Δ = max|y−ȳ| on T. With
p = 5 the vector radius is ~10⁻⁷ against eigenvalue margins ~10⁻⁴.

**Matrix check.** N_h − D_h ⪰ 0 for *all* interval data reduces, by Lemma R
and ‖ΔQ‖₂ ≤ max row sum of the radius matrix, to λ_min(Â) > ε_tot with Â
the float midpoint combination (safe coefficient directions) and ε_tot the
sum of: matrix-radius norms, rank-one radius shifts, and the float build
error 8u·max row sum of the absolute-term matrix.

λ_min(Â) > ε_tot is certified by the method of **Rump, *Verification of
positive definiteness*, BIT Numer. Math. 46 (2006) 433–452**, implemented
literally (`rump_psd_certificate` in `m3_certify.py`):

- *Theorem 2.3* (op. cit.): if the floating-point Cholesky factorization of
  a symmetric A ∈ 𝔽ⁿˣⁿ runs to completion, then λ_min(A) > −‖Δ(A)‖₂, where
  Δ(A)_ij := α_ij d_i d_j + M·eta, α_ij := γ_{s(i,j)+2} ≤ γ_{n+1},
  d_j := ((1−α_jj)⁻¹ a_jj)^{1/2}, M := 3(2n + max a_νν), with
  eps = 2⁻⁵³, eta = 2⁻¹⁰⁷⁴, γ_k = k·eps/(1−k·eps) (IEEE-754 double,
  round-to-nearest; underflow allowed). We bound
  ‖Δ(A)‖₂ ≤ γ_{n+1}·Σ_j d_j² + n·M·eta (entrywise domination + Perron).
- *Corollary 2.4 + Lemma 2.5* (op. cit.): with c ≥ ‖Δ(A)‖₂ + ε_tot in 𝔽,
  set ã_ii := fl(dᵢ′ − ϕ|dᵢ′|), dᵢ′ := fl(a_ii − c), ϕ := eps(1+2eps)
  (this guarantees ã_ii ≤ a_ii − c without directed rounding). If floating
  Cholesky of Ã runs to completion, the corollary's proof chain gives
  λ_min(A) ≥ c + λ_min(Ã) > c − ‖Δ(Ã)‖₂ ≥ ε_tot, as required.
- Library routine: Rump, op. cit., p. 1–2: "A major advantage of the method
  is that any library routine can be used." The code calls LAPACK dpotrf
  (scipy); this is the same usage as INTLAB's `isspd`.
- The float evaluation of the scalar bound c itself is padded by explicit
  (1 + n·eps)-type inflation factors, recorded in code comments.

**Frozen certificate parameters** (m3_certify.py defaults): mesh 12×6×6
(2592 tets, 5544 CR dofs), Y = 1.25, ρ = 55, θ = 0.7, θ₂ = θ' = α = 0.5,
ν* = 1.05, 8 uniform s-windows, Taylor p = 5, arb prec 128. Result
(2026-07-10): all windows certified.

## 6. Residual open items (citation-level)

- ~~κ₁ hypotheses~~ **resolved 2026-07-10**: the certificate no longer
  depends on CZZ20 — Lemma I1 (self-contained, arbitrary tets,
  κ_sc = √(1/π²+1/15), inputs: Payne–Weinberger 1960/Bebendorf 2003) is the
  default; CZZ20's κ₁ = √(1/π²+1/120) is an optional sharpening
  (`kappa_mode="czz"`), to be audited only if ever needed.
- ~~floating-Cholesky constant~~ **resolved 2026-07-10**: exact
  implementation of Rump BIT 46 (2006) Thm 2.3 + Cor 2.4 + Lemma 2.5
  (see §5); the earlier 10× guard is retired.
- [TODO] CAMWA 2026 paper into `refs/`; reconcile/cite as primary LEB source.
- G4 of DESIGN.md: EGM Ch. 6 (continuous spectrum [1,∞), discreteness
  below), Shimizu (Lemma 0), Aronszajn/real-analyticity (unique
  continuation).
