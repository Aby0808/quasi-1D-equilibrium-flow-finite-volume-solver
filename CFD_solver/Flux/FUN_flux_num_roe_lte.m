function [F, meta] = FUN_flux_num_roe_lte(UL, UR)
% Roe flux-difference splitting for real-gas LTE, based on:
% "A simple extension of Roe’s scheme for real gases" (Arabi et al., JCP 2017).
% Adds an entropy fix on |lambda| to avoid expansion shocks.

% Ensure column vectors
UL = UL(:);
UR = UR(:);

% ----------------------------
% 1) Left/right primitive + thermo
% ----------------------------
[rhoL, uL, PL, HL, cL] = local_state(UL);
[rhoR, uR, PR, HR, cR] = local_state(UR);

% Physical fluxes
FL = [ UL(2);
       UL(2)^2/UL(1) + PL;
       (UL(2)/UL(1))*(UL(3) + PL) ];

FR = [ UR(2);
       UR(2)^2/UR(1) + PR;
       (UR(2)/UR(1))*(UR(3) + PR) ];

dF = FR - FL;

% ----------------------------
% 2) Roe averages (same operator as ideal gas)
% ----------------------------
sL = sqrt(rhoL);
sR = sqrt(rhoR);
den = sL + sR;

rho_t = sL * sR;  % = sqrt(rhoL*rhoR)
u_t   = (sL*uL + sR*uR) / den;
H_t   = (sL*HL + sR*HR) / den;

% Sound speed choice (Method III in the paper: Roe-average of c)
c_t   = (sL*cL + sR*cR) / den;
c_t   = max(c_t, 1e-12);  % safety

% ----------------------------
% 3) Eigenvalues and wave strengths (alphas)
% ----------------------------
dp   = PR - PL;
du   = uR - uL;
drho = rhoR - rhoL;

lam1 = u_t + c_t;
lam2 = u_t - c_t;
lam3 = u_t;

ct2  = c_t^2;

alpha1 = (dp + rho_t*c_t*du) / (2*ct2);
alpha2 = (dp - rho_t*c_t*du) / (2*ct2);
alpha3 = drho - dp/ct2;

% Eigenvectors (only need e1,e2 fully; e3 third component is avoided by the paper trick)
e1 = [1; u_t + c_t; H_t + u_t*c_t];
e2 = [1; u_t - c_t; H_t - u_t*c_t];

% ----------------------------
% 3.5) Entropy fix: smooth |lambda| near sonic points
% ----------------------------
% Typical choice: delta = epsEntropy * c_t (tune epsEntropy in [0.05, 0.2])
epsEntropy = 0.10;
delta = max(epsEntropy * c_t, 1e-12);

absl1 = abs_entropy(lam1, delta);
absl2 = abs_entropy(lam2, delta);
absl3 = abs_entropy(lam3, delta);

% ----------------------------
% 4) Compute |ΔF| with Arabi et al. non-singular trick
% ----------------------------
absdF = zeros(3,1);

% Component 1
absdF(1) = absl1*alpha1 + absl2*alpha2 + absl3*alpha3;

% Component 2
absdF(2) = (u_t + c_t)*absl1*alpha1 + (u_t - c_t)*absl2*alpha2 + u_t*absl3*alpha3;

% Component 3: special treatment (Eq. 44–46 in the paper)
X    = dF(3) - e1(3)*lam1*alpha1 - e2(3)*lam2*alpha2;
absX = X * sign_nonzero(u_t);

absdF(3) = e1(3)*absl1*alpha1 + e2(3)*absl2*alpha2 + absX;

% ----------------------------
% 5) Roe numerical flux
% ----------------------------
F = 0.5*(FL + FR) - 0.5*absdF;

% Diagnostics
meta.rhoL = rhoL; meta.uL = uL; meta.PL = PL; meta.HL = HL; meta.cL = cL;
meta.rhoR = rhoR; meta.uR = uR; meta.PR = PR; meta.HR = HR; meta.cR = cR;
meta.rho_t = rho_t; meta.u_t = u_t; meta.H_t = H_t; meta.c_t = c_t;
meta.lam = [lam1; lam2; lam3];
meta.alpha = [alpha1; alpha2; alpha3];
meta.absLambda = [absl1; absl2; absl3];
meta.entropy.delta = delta;
meta.entropy.epsEntropy = epsEntropy;
meta.absdF = absdF;

end

% ============================================================
% Helpers
% ============================================================

function [rho, u, P, H, c] = local_state(U)
rho = U(1);
m   = U(2);
Et  = U(3);

%fprintf('\n%d%d%d',U(1),U(2),U(3))

u = m / rho;

E = Et / rho;
e = E - 0.5*u*u;  % internal energy per mass

st = FUN_get_eq_prop('from_rhoe', rho, e);

P = st.P;

% total enthalpy per mass
H = (Et + P) / rho;

% sound speed (starting point): c = sqrt(gamma*P/rho)
c = sqrt(max(st.gamma * P / rho, 0));
end

function s = sign_nonzero(x)
% Prevents returning 0 at exactly x=0 (important for |X| = X*sign(u_t)).
if x >= 0
    s = 1;
else
    s = -1;
end
end

function a = abs_entropy(lam, delta)
% Harten-type entropy fix for Roe:
% |lam|_e = |lam|,                    if |lam| >= delta
%        = 0.5*(lam^2/delta + delta), if |lam| <  delta
x = abs(lam);
if x >= delta
    a = x;
else
    a = 0.5*(x*x/delta + delta);
end
end
