import Mathlib
open MeasureTheory ProbabilityTheory Real Set intervalIntegral

noncomputable section
namespace FbmCore

/-- Left-sided Riemann-Liouville integral with lower terminal 0. -/
def RL (α : ℝ) (f : ℝ → ℝ) (t : ℝ) : ℝ :=
  (Gamma α)⁻¹ * ∫ s in (0:ℝ)..t, (t - s) ^ (α - 1) * f s

/-- Variance `u^{2H}` as an `NNReal`. -/
def var (H u : ℝ) : NNReal := ⟨|u| ^ (2 * H), rpow_nonneg (abs_nonneg u) _⟩

@[simp] lemma var_coe (H u : ℝ) : (var H u : ℝ) = |u| ^ (2 * H) := rfl

lemma gaussian_at_zero (H u : ℝ) (hu : 0 < u) :
    gaussianPDFReal 0 (var H u) 0 = (√(2 * π))⁻¹ * u ^ (-H) := by
  have hs : √(u ^ (2 * H)) = u ^ H := by
    rw [show 2 * H = H * 2 by ring, rpow_mul hu.le, rpow_two]
    exact Real.sqrt_sq (rpow_nonneg hu.le _)
  rw [gaussianPDFReal_def]
  simp only [var_coe, sub_zero, abs_of_pos hu]
  rw [show (0:ℝ) ^ 2 = 0 by norm_num, neg_zero, zero_div, Real.exp_zero, mul_one,
    Real.sqrt_mul (by positivity), hs, rpow_neg hu.le, mul_inv]

lemma kernel_reversal (f : ℝ → ℝ) (α t : ℝ) :
    ∫ u in (0:ℝ)..t, f (t - u) * u ^ (α - 1) = ∫ s in (0:ℝ)..t, (t - s) ^ (α - 1) * f s := by
  have h := intervalIntegral.integral_comp_sub_left (fun s => (t - s) ^ (α - 1) * f s) (a := 0) (b := t) t
  simp only [sub_sub_cancel, sub_self, sub_zero] at h
  rw [← h]
  congr 1; ext u; ring

lemma integral_power (α t : ℝ) (hα : 0 < α) :
    ∫ u in (0:ℝ)..t, u ^ (α - 1) = t ^ α / α := by
  rw [integral_rpow (Or.inl (by linarith))]
  simp [hα.ne']

lemma normalized_kernel_mass (α t : ℝ) (hα : 0 < α) (ht : 0 < t) :
    ∫ s in (0:ℝ)..t, α * t ^ (-α) * (t - s) ^ (α - 1) = 1 := by
  have h := kernel_reversal (fun _ => (1:ℝ)) α t
  simp only [one_mul, mul_one] at h
  rw [intervalIntegral.integral_const_mul, ← h, integral_power α t hα, rpow_neg ht.le]
  have : t ^ α ≠ 0 := (rpow_pos_of_pos ht α).ne'
  field_simp


/-- Step 2 of the proof (algebraic core): adding nonnegative diagonal noise to a
2x2 covariance matrix can only increase its determinant. -/
lemma det2_mono (a b c e d : ℝ) (ha : 0 ≤ a) (hc : 0 ≤ c) (he : 0 ≤ e) (hd : 0 ≤ d) :
    a * c - b ^ 2 ≤ (a + e) * (c + d) - b ^ 2 := by
  nlinarith [mul_nonneg ha hd, mul_nonneg he hc, mul_nonneg he hd]

/-- The expected total mass: `∫_0^t (2π)^{-1/2} u^{α-1} du = t^α / (α √(2π))`. -/
lemma expected_mass (α t : ℝ) (hα : 0 < α) :
    ∫ u in (0:ℝ)..t, (√(2 * π))⁻¹ * u ^ (α - 1) = t ^ α / (α * √(2 * π)) := by
  rw [intervalIntegral.integral_const_mul, integral_power α t hα]
  field_simp

/-! ## The probabilistic part, with the fBm input as a labeled hypothesis -/

open scoped ENNReal

/-- Gaussian density at `0` for variance `|u|^{2H}`, as an extended nonnegative real. -/
def dens (H u : ℝ) : ℝ≥0∞ := ENNReal.ofReal (gaussianPDFReal 0 (var H u) 0)

lemma measurable_dens (H : ℝ) : Measurable (dens H) := by
  unfold dens
  apply Measurable.ennreal_ofReal
  simp only [gaussianPDFReal_def, var_coe]
  fun_prop

lemma dens_eq {α H u : ℝ} (hH : H = 1 - α) (hu : 0 < u) :
    dens H u = ENNReal.ofReal ((√(2 * π))⁻¹ * u ^ (α - 1)) := by
  rw [dens, gaussian_at_zero H u hu, hH, show -(1 - α) = α - 1 by ring]

section Prob
variable {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω) (L : Kernel Ω ℝ)

/-- **Hypothesis (H-int).** The intensity of the random measure `L` (the local time at `0`
of the fBm on `[0,t]`, in the forward time variable `u`) has density `p_u(0)`, the Gaussian
density at `0` with variance `u^{2H}`. In the paper this is equation (4); it follows from the
occupation density formula, the Gaussian law of `B_u`, and uniform integrability
(Steps 1-3). Mathlib has no fBm or local time, so it is assumed here. -/
def IntensityHyp (H t : ℝ) : Prop :=
  ∀ A : Set ℝ, MeasurableSet A → ∫⁻ ω, L ω A ∂P = ∫⁻ u in A ∩ Ioc 0 t, dens H u

lemma intensity_eq {H t : ℝ} (h : IntensityHyp P L H t) :
    L ∘ₘ P = (volume.restrict (Ioc 0 t)).withDensity (dens H) := by
  ext A hA
  rw [Measure.bind_apply hA L.measurable.aemeasurable, h A hA, withDensity_apply _ hA,
    Measure.restrict_restrict hA]

/-- Step 4 for nonnegative functions, without time reversal. -/
theorem expected_lintegral {α H t : ℝ} (hH : H = 1 - α) (h : IntensityHyp P L H t)
    (g : ℝ → ℝ≥0∞) (hg : Measurable g) :
    ∫⁻ ω, ∫⁻ u, g u ∂(L ω) ∂P
      = ∫⁻ u in Ioc 0 t, ENNReal.ofReal ((√(2 * π))⁻¹ * u ^ (α - 1)) * g u := by
  rw [← Measure.lintegral_bind L.measurable.aemeasurable hg.aemeasurable, intensity_eq P L h,
    lintegral_withDensity_eq_lintegral_mul _ (measurable_dens H) hg]
  apply setLIntegral_congr_fun measurableSet_Ioc
  intro u hu
  simp only [Pi.mul_apply, dens_eq hH hu.1]

/-- Step 4 with time reversal `s = t - u`, nonnegative functions. -/
theorem expected_lintegral_reversed {α H t : ℝ} (hH : H = 1 - α) (h : IntensityHyp P L H t)
    (f : ℝ → ℝ≥0∞) (hf : Measurable f) :
    ∫⁻ ω, ∫⁻ u, f (t - u) ∂(L ω) ∂P
      = ∫⁻ s in Ioc 0 t, ENNReal.ofReal ((√(2 * π))⁻¹ * (t - s) ^ (α - 1)) * f s := by
  rw [expected_lintegral P L hH h (fun u => f (t - u)) (hf.comp (measurable_const.sub measurable_id))]
  set F : ℝ → ℝ≥0∞ := fun u => ENNReal.ofReal ((√(2 * π))⁻¹ * u ^ (α - 1)) * f (t - u) with hF
  have hFm : Measurable F := by
    simp only [hF]
    apply Measurable.mul
    · apply Measurable.ennreal_ofReal; fun_prop
    · exact hf.comp (measurable_const.sub measurable_id)
  rw [← lintegral_indicator measurableSet_Ioc,
    ← lintegral_sub_left_eq_self (fun u => (Ioc 0 t).indicator F u) t]
  have key : (fun s => (Ioc 0 t).indicator F (t - s)) = (Ico 0 t).indicator (fun s => F (t - s)) := by
    ext s
    by_cases hs : s ∈ Ico 0 t
    · have : t - s ∈ Ioc 0 t := ⟨by linarith [hs.2], by linarith [hs.1]⟩
      simp [indicator_of_mem hs, indicator_of_mem this]
    · have : t - s ∉ Ioc 0 t := by
        intro h'; apply hs; exact ⟨by linarith [h'.2], by linarith [h'.1]⟩
      simp [indicator_of_notMem hs, indicator_of_notMem this]
  rw [key, lintegral_indicator measurableSet_Ico]
  rw [setLIntegral_congr Ico_ae_eq_Ioc]
  apply setLIntegral_congr_fun measurableSet_Ioc
  intro s _
  simp [hF, sub_sub_cancel]

/-- **Main theorem (formal version).** Under hypothesis (H-int), for every Borel `f` with
`∫_0^t |f(s)| (t-s)^{α-1} ds < ∞`, the integral of `f(t-u)` against `L` is integrable
for almost every path and
`E ∫ f(t-u) L(du) = Γ(α)/√(2π) · I^α f(t)`. -/
theorem main_theorem [SFinite P] [IsSFiniteKernel L] {α H t : ℝ} (hα : 0 < α) (hH : H = 1 - α)
    (ht : 0 ≤ t) (h : IntensityHyp P L H t) (f : ℝ → ℝ) (hf : Measurable f)
    (hC : IntegrableOn (fun s => (t - s) ^ (α - 1) * f s) (Ioc 0 t)) :
    (∀ᵐ ω ∂P, Integrable (fun u => f (t - u)) (L ω)) ∧
    ∫ ω, ∫ u, f (t - u) ∂(L ω) ∂P = Gamma α / √(2 * π) * RL α f t := by
  set F : ℝ → ℝ := fun u => f (t - u) with hF
  have hFm : Measurable F := hf.comp (measurable_const.sub measurable_id)
  -- integrability with respect to the intensity measure
  have hInt : Integrable F (L ∘ₘ P) := by
    refine ⟨hFm.aestronglyMeasurable, ?_⟩
    rw [hasFiniteIntegral_iff_eNorm,
      Measure.lintegral_bind L.measurable.aemeasurable hFm.eNorm.aemeasurable]
    have := expected_lintegral_reversed P L hH h (fun s => ‖f s‖ₑ) hf.eNorm
    simp only [hF] at this ⊢
    rw [this]
    have hfin := (hC.const_mul (√(2 * π))⁻¹).hasFiniteIntegral
    rw [hasFiniteIntegral_iff_eNorm] at hfin
    refine lt_of_le_of_lt (le_of_eq ?_) hfin
    apply setLIntegral_congr_fun measurableSet_Ioc
    intro s hs
    have h1 : 0 ≤ (√(2 * π))⁻¹ * (t - s) ^ (α - 1) :=
      mul_nonneg (inv_nonneg.2 (Real.sqrt_nonneg _)) (rpow_nonneg (by linarith [hs.2]) _)
    dsimp only
    rw [← mul_assoc, eNorm_mul, Real.eNorm_of_nonneg h1]
  refine ⟨Measure.ae_integrable_of_integrable_comp hInt, ?_⟩
  have step1 : ∫ ω, ∫ u, F u ∂(L ω) ∂P = ∫ u, F u ∂(L ∘ₘ P) := by
    rw [Measure.comp_eq_comp_const_apply] at hInt ⊢
    rw [Kernel.integral_comp hInt]
    simp [Kernel.const_apply]
  rw [step1, intensity_eq P L h,
    integral_withDensity_eq_integral_toReal_smul (measurable_dens H)
      (ae_of_all _ (fun _ => ENNReal.ofReal_lt_top))]
  have step2 : ∫ u in Ioc 0 t, (dens H u).toReal • F u
      = ∫ u in Ioc 0 t, (√(2 * π))⁻¹ * (F u * u ^ (α - 1)) := by
    apply setIntegral_congr_fun measurableSet_Ioc
    intro u hu
    have h1 : 0 ≤ (√(2 * π))⁻¹ * u ^ (α - 1) :=
      mul_nonneg (inv_nonneg.2 (Real.sqrt_nonneg _)) (rpow_nonneg hu.1.le _)
    simp only [dens_eq hH hu.1, ENNReal.toReal_ofReal h1, smul_eq_mul]
    ring
  rw [step2, MeasureTheory.integral_const_mul, ← intervalIntegral.integral_of_le ht]
  simp only [hF]
  rw [kernel_reversal f α t, RL]
  have hG : Gamma α ≠ 0 := (Gamma_pos_of_pos hα).ne'
  field_simp

end Prob

end FbmCore

