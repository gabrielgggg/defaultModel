# Sovereign Default Model with Long-Term Debt
Computed with taste shocks (discrete choice methods) for $B'$ and $d$. Fortran parallelized with OpenMP.

Value and default option:

$$ V\left(y, B\right) = \rho_D \log \left\\{ \exp\left[ \dfrac{V^d(y)}{\rho_D} \right] + \exp \left[ \dfrac{V^r(y, B)}{\rho_D} \right] \right\\} $$

$$ \Pr\left(d = 1 \middle| y, B\right) = \dfrac{   \exp\left[ \dfrac{V^d(y)}{\rho_D} \right] }{\exp\left[ \dfrac{V^d(y)}{\rho_D} \right] + \exp \left[ \dfrac{V^r(y, B)}{\rho_D} \right]} = \dfrac{1}{1 + \exp \left[ \dfrac{V^r(y, B) - V^d(y)}{\rho_D} \right]} $$

Default value:

$$ V^d\left(y\right) = u\left[h\left(y\right)\right] + \beta \mathbb{E}_{y'|y} \left\\{ \gamma V\left( y', 0 \right) + (1-\gamma) V^d\left(y'\right) \right\\} $$

Repayment values:

$$ W\left(y, B, B'\right) = u\left[ y - \kappa B + q\left(y, B'\right) \left( B' - (1-\delta) B \right) \right] + \beta \mathbb{E}_{y'|y} V\left(y', B'\right) $$

$$ V^r\left(y, B\right) = \rho_B \log \sum_{B'} \exp \left[ \dfrac{W\left(y, B, B'\right))}{\rho_B} \right] $$

Choice probabilities for $B'$:

$$ \Pr\left(B' = x \middle| y, B \right) = \dfrac{ \exp \left[ \dfrac{W\left(y, B, x\right))}{\rho_B} \right] }{\sum_{i} \exp \left[ \dfrac{W\left(y, B, i\right))}{\rho_B} \right]} $$

Bond price schedule:

$$ q\left(y, B'\right) = \dfrac{1}{1+r} \mathbb{E}_{y'|y} \Pr\left(d=0 \middle| y', B' \right) \left[ \kappa + (1-\delta) \mathcal{Q}(y', B') \right] $$

with

$$ \mathcal{Q}(y', B') = \sum_{B''} \Pr\left( B'' \middle| y', B' \right) q\left(y', B''\right) $$

Functional forms and shocks:

$$ u\left(c\right) = \dfrac{c^{1-\sigma} - 1}{1-\sigma} $$

$$ h\left(y\right) = y - \max\\{ 0, \lambda_0 y + \lambda_1 y^2 \\} $$

$$ \log y' = - (1-\rho) \dfrac{\sigma_y^2}{2 \left( 1 - \rho^2\right)} + \rho \log y + \sigma_y \varepsilon, \quad \varepsilon \sim \mathcal{N}(0, 1) $$
