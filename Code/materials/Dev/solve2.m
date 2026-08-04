function [p, q, peqp, sigma, stif, flag, iter, res] = ...
    solve2(sigma0, devp, peqp0, ni, lambda, kappa, mu, tau, cohes, M, fric_ang, dt, itmax, atol, rtol)

% Identity tensor in vector form
m = [1;1;1;0;0;0];
% Scaling matrix for strain in Voigt notation
Ds = diag([1;1;1;0.5;0.5;0.5]);
% Deviatoric operator in matrix form
Dev = eye(6) - (1/3)*(m*m');

p0 = max((sigma0(1)+sigma0(2)+sigma0(3))/3, 1e-10);
s0 = sigma0;
s0(1:3) = s0(1:3) - p0;

Kbar = 1 / kappa;
Gbar = (3/2) * (1 - 2*ni) / (1 + ni) * Kbar;

devpv = devp(1) + devp(2) + devp(3);
sbar = zeros(6,1);
sbar(1:3) = 2 * Gbar * (devp(1:3) - devpv/3);
sbar(4:6) = Gbar * devp(4:6);
pbar = Kbar * (m'*devp);

sbar_sbar = sbar(1)^2 + sbar(2)^2 + sbar(3)^2 + 2*(sbar(4)^2 + sbar(5)^2 + sbar(6)^2);
sbar_s0   = sbar(1)*s0(1) + sbar(2)*s0(2) + sbar(3)*s0(3) + 2*(sbar(4)*s0(4) + sbar(5)*s0(5) + sbar(6)*s0(6));
s0_s0     = s0(1)^2 + s0(2)^2 + s0(3)^2 + 2*(s0(4)^2 + s0(5)^2 + s0(6)^2);

% Initialize solution values (p,q) - Guarded against negative round-off
p = p0 * pbar + p0;
if p < 1e-10, p = p0; end
q = sqrt(1.5 * max(p0^2*sbar_sbar + 2*p0*sbar_s0 + s0_s0, 0.0));

iter = 0;
conv = false;
condIssue = false;
res = zeros(itmax,1);

% Convert degrees to radians for tan()
tanPhi = tan(deg2rad(fric_ang));
shift = cohes / max(tanPhi, 1e-12);

muSafe = max(mu, 1e-12);

while (~conv)
    iter = iter + 1;

    phat = max(p + shift, 1e-10);
    peq = phat + q^2 / (M^2 * phat);
    dpeq_dp = 1 - q^2 / (M^2 * phat^2);
    dpeq_dq = 2*q / (M^2 * phat);
    d2peq_dp2 = 2*q^2 / (M^2 * phat^3);
    d2peq_dpdq = -2*q / (M^2 * phat^2);
    d2peq_dq2 = 2 / (M^2 * phat);

    ptrial = pbar*p + p0;
    qtrial = sqrt(1.5 * max(p^2*sbar_sbar + 2*p*sbar_s0 + s0_s0, 0.0));

    K = Kbar * p;
    G = Gbar * p;
    
    qtrialSafe = max(qtrial, 1e-10);
    Devp = (qtrial - q) / (3*G);
    Devpv = (ptrial - p) / K;

    dDevp_dp = 1 / (2*G*qtrialSafe) * (p*sbar_sbar + sbar_s0) - Devp/p;
    dDevp_dq = -1.0 / (3*G);
    dDevpv_dp = -p0 / (K*p);

    peqp = peqp0 * exp(Devpv / (lambda - kappa));
    dpeqp_dp = peqp * dDevpv_dp / (lambda - kappa);

    psi1 = dpeq_dp * Devp - dpeq_dq * Devpv;
    
    peqpSafe = max(peqp, 1e-12);
    ratio = min(max(peq / peqpSafe, 1e-12), 1e10);
    exponent = (lambda - kappa) / muSafe;
    psi2 = Devpv - (mu / tau) * dt * (ratio^exponent);

    dpsi1_dp = d2peq_dp2*Devp + dpeq_dp*dDevp_dp - d2peq_dpdq*Devpv - dpeq_dq*dDevpv_dp;
    dpsi1_dq = d2peq_dpdq*Devp + dpeq_dp*dDevp_dq - d2peq_dq2*Devpv;

    fac = ((lambda - kappa) / tau) * dt * (ratio^(exponent - 1)) / peqpSafe;
    dpsi2_dp = dDevpv_dp - fac*dpeq_dp + fac*(peq/peqpSafe)*dpeqp_dp;
    dpsi2_dq = -fac*dpeq_dq;

    J = [dpsi1_dp, dpsi1_dq;
         dpsi2_dp, dpsi2_dq];

    r = [psi1; psi2];

    colNorms = sqrt(max(diag(J'*J), 1e-20));
    Dc = diag(1 ./ colNorms);
    Jscaled = J * Dc;

    if (rcond(Jscaled) < 1e-15)
        condIssue = true;
    else
        dx = -Jscaled \ r;
        dx = Dc * dx;

        p = max(p + dx(1), 1e-10);
        q = max(q + dx(2), 0.0);
    end

    rnorm = norm(r);
    if (iter == 1)
        rnorm0 = max(rnorm, 1.0);
    end
    res(iter) = rnorm;

    conv = ((rnorm < rtol*rnorm0 + atol) || iter >= itmax) || condIssue;
end

flag = (rnorm < rtol*rnorm0 + atol) && ~condIssue;
res = res(1:iter);

% Direction tensor update (zero vector when qtrial -> 0)
if (qtrial > 1e-12)
    strial = p*sbar + s0;
    nt = 1.5 * strial / qtrial;
else
    nt = zeros(6,1);
end

% Reconstructed stress tensor
sigma = (2/3)*q*nt + p*m;

% Algorithmic Tangent Calculation
phat = max(p + shift, 1e-10);
peq = phat + q^2 / (M^2 * phat);
dpeq_dp = 1 - q^2 / (M^2 * phat^2);
dpeq_dq = 2*q / (M^2 * phat);
d2peq_dp2 = 2*q^2 / (M^2 * phat^3);
d2peq_dpdq = -2*q / (M^2 * phat^2);
d2peq_dq2 = 2 / (M^2 * phat);

dptrial_dp = pbar;
qtrialSafe = max(qtrial, 1e-10);
dqtrial_dp = 1.5 / qtrialSafe * (p*sbar_sbar + sbar_s0);

peqpSafe = max(peqp, 1e-12);
B = dt * ((lambda - kappa) / tau) * (min(max(peq / peqpSafe, 1e-12), 1e10)^((lambda - kappa)/muSafe - 1)) / peqpSafe;
A = 1 + B * (peq / peqpSafe) * (peqp0 / (lambda - kappa)) * exp(Devpv / (lambda - kappa));

% Prevent division by zero at Critical State Line (dpeq_dp = 0)
if abs(dpeq_dp) < 1e-4
    betaSafe = 1e-4;
else
    betaSafe = dpeq_dp;
end
dlambda = Devpv / betaSafe;

F11 = 1 + K*dpeq_dp*(B/A) + Kbar*dpeq_dp*dlambda - dptrial_dp;
F12 = K*dpeq_dq*(B/A);
F21 = 3*G*d2peq_dpdq*dlambda + 3*G*dpeq_dq*(B/A) - 3*G*(dpeq_dq/betaSafe)*d2peq_dp2*dlambda + 3*Gbar*dpeq_dq*dlambda - dqtrial_dp;
F22 = 1 + 3*G*d2peq_dq2*dlambda + 3*G*(dpeq_dq^2/betaSafe)*(B/A) - 3*G*(dpeq_dq/betaSafe)*d2peq_dpdq*dlambda;

detF = F11*F22 - F12*F21;
if abs(detF) < 1e-14
    if detF < 0
        detF = -1e-14;
    else
        detF = 1e-14;
    end
end

c11 =  F22 / detF;
c12 = -F12 / detF;
c21 = -F21 / detF;
c22 =  F11 / detF;

edev = nt' * devp;
t = sbar / (2*Gbar);

d1 = (2/3) * (c21 - 2*Gbar*(q/qtrialSafe)*c11*edev) * K;
d2 = 2 * G * c12;
d3 = (4/3) * ((c22 - q/qtrialSafe) - 2*Gbar*(q/qtrialSafe)*c12*edev) * G;
d4 = c11 * K;
d5 = 2 * G * (q/qtrialSafe);
d6 = 2 * Gbar * (q/qtrialSafe) * c11 * K;
d7 = 4 * Gbar * (q/qtrialSafe) * c12 * G;

stif = d1*(nt*m') + d2*(m*nt') + d3*(nt*nt') + d4*(m*m') + d5*Dev*Ds + d6*(t*m') + d7*(t*nt');

% Symmetrize tangent matrix for PCG linear solver compatibility
stif = 0.5 * (stif + stif');

end