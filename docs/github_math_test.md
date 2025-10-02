# GitHub Math Rendering Test

This file demonstrates how different LaTeX math formats appear on GitHub.

## 1. Inline Math (Single Dollar Signs) - ❓ Support Varies

**What you write:**
```
Variables like $\alpha$, $\beta$, and $T_{surf}$ should render as math.
Units like m$^{2}$s$^{-2}$ should show superscripts.
```

**How it appears on GitHub:**
Variables like $\alpha$, $\beta$, and $T_{surf}$ should render as math.
Units like m$^{2}$s$^{-2}$ should show superscripts.

## 2. Block Math (Double Dollar Signs) - ✅ Supported Since 2022

**What you write:**
```
$$
\langle \bar{\varphi \;} \rangle \left(z\right)=\frac{1}{A_f }\int_{{\;}_{\Omega_{f\;} } } \varphi \;\textrm{dA}
$$
```

**How it appears on GitHub:**
$$
\langle \bar{\varphi \;} \rangle \left(z\right)=\frac{1}{A_f }\int_{{\;}_{\Omega_{f\;} } } \varphi \;\textrm{dA}
$$

## 3. Examples from Your Tutorial Files

**These inline math expressions from your auto-generated tutorials:**
```matlab
format_surface_plot('$\alpha$ [-]')
xlabel('$x$ [m]','Interpreter','latex')  
ylabel('$T$ [K]','Interpreter','latex')
format_surface_plot('$p$ [m$^{2}$s$^{-2}$]')
format_surface_plot('$K^*$ [Wm$^{-2}$]')
leg{n} = ['$t=', num2str(tsel(n), '%8.0f'), '$ s'];
```

**Will likely show as literal text:**
- `$\alpha$` instead of α  
- `m$^{2}$s$^{-2}$` instead of m²s⁻²
- `$K^*$` instead of K*

## 4. GitHub-Compatible Alternatives

**✅ Unicode symbols (guaranteed to work):**
- α β γ δ θ λ μ ρ σ τ φ ω (Greek letters)  
- ⁰ ¹ ² ³ ⁴ ⁵ ⁶ ⁷ ⁸ ⁹ (superscripts)
- ₀ ₁ ₂ ₃ ₄ ₅ ₆ ₇ ₈ ₉ (subscripts)
- ⁻ (minus for negative powers)

**Examples of conversions:**
- `$\alpha$` → `α`  
- `m$^{2}$s$^{-2}$` → `m²s⁻²`
- `$T_{surf}$` → `Tₛᵤᵣf` or `T_surf`
- `$K^*$` → `K*` or `K⁎`
- `Wm$^{-2}$` → `Wm⁻²`

## 5. Current Status in Your Repository

Check these files on GitHub to see actual rendering:
- [docs/udales-fields-tutorial.md](https://github.com/uDALES/u-dales/blob/enhance-udbase-udgeom-integration/docs/udales-fields-tutorial.md#L160) (has block math)
- Your auto-generated tutorial files (will have inline math issues)

**Expected issues:** Most `$...$` expressions in MATLAB code blocks will show as literal dollar signs rather than rendered math.