function [p, q, peqp, sigma, stif, flag, iter, res] = ...
    solve2(sigma0, devp, peqp0, ni, lambda, kappa, mu, tau, cohes, M, fric_ang, dt, itmax, atol, rtol)

% Identity tensor in vector form
m = [1;1;1;0;0;0];
% Scaling matrix for strain in Voight notation
Ds = diag([1;1;1;0.5;0.5;0.5]);
% Deviatoric operator in matrix form
Dev = eye(6) - 1/3*(m*m');

p0 = (sigma0(1)+sigma0(2)+sigma0(3))/3;
% s0 = Dev*sigma0;
s0 = sigma0;
s0(1:3) = s0(1:3) - p0;

Kbar = 1/kappa;
Gbar = 3/2*(1-2*ni)/(1+ni)*Kbar;
% sbar = 2*Gbar*(Dev*(Ds*devp));
devpv = devp(1)+devp(2)+devp(3);
sbar = zeros(6,1);
sbar(1:3) = 2*Gbar*(devp(1:3) - devpv/3);
sbar(4:6) = Gbar*devp(4:6);
pbar = Kbar*(m'*devp);
sbar_sbar = sbar(1)*sbar(1) + sbar(2)*sbar(2) + sbar(3)*sbar(3) + 2*(sbar(4)*sbar(4) + sbar(5)*sbar(5) + sbar(6)*sbar(6));
sbar_s0 = sbar(1)*s0(1) + sbar(2)*s0(2) + sbar(3)*s0(3) + 2*(sbar(4)*s0(4) + sbar(5)*s0(5) + sbar(6)*s0(6));
s0_s0 = s0(1)*s0(1) + s0(2)*s0(2) + s0(3)*s0(3) + 2*(s0(4)*s0(4) + s0(5)*s0(5) + s0(6)*s0(6));

% Initialize solution values (p,q)
p = p0*pbar + p0;
q = sqrt(3/2*(p0^2*sbar_sbar + 2*p0*sbar_s0 + s0_s0));

iter = 0;
conv = false;
condIssue = false;
res = zeros(itmax,1);
while (~conv)

    iter = iter + 1;

    phat = p + cohes/tan(fric_ang);
    peq = phat + q^2/(M^2*phat);
    dpeq_dp = 1 - q^2/(M^2*phat^2);
    dpeq_dq = 2*q/(M^2*phat);
    d2peq_dp2 = 2*q^2/(M^2*phat^3);
    d2peq_dpdq = -2*q/(M^2*phat^2);
    d2peq_dq2 = 2/(M^2*phat);

    ptrial = pbar*p + p0;
    qtrial = sqrt(3/2*(p^2*sbar_sbar + 2*p*sbar_s0 + s0_s0));

    K = Kbar*p;
    G = Gbar*p;
    Devp = (qtrial - q)/(3*G);
    Devpv = (ptrial - p)/K;

    dDevp_dp = 1/(2*G*qtrial)*(p*sbar_sbar + sbar_s0) - Devp/p;
    dDevp_dq = -1.0/(3*G);
    dDevpv_dp = -p0/(K*p);

    peqp = peqp0*exp(Devpv/(lambda-kappa));
    %peqp = peqp0/(1.0-Devpv/(lambda-kappa));
    dpeqp_dp = peqp0*exp(Devpv/(lambda-kappa))*dDevpv_dp/(lambda-kappa);

    % Using R1 = dpeq_dp*Devp - dpeq_dq*Devpv
    % instead of R1 = Devp - dpeq_dq/dpeq_dp*Devpv
    % because dpeq_dp can be zero
    psi1 = dpeq_dp*Devp - dpeq_dq*Devpv;
    psi2 = Devpv - mu/tau*dt*(peq/peqp)^((lambda-kappa)/mu);

    dpsi1_dp = d2peq_dp2*Devp + dpeq_dp*dDevp_dp - d2peq_dpdq*Devpv - dpeq_dq*dDevpv_dp;
    dpsi1_dq = d2peq_dpdq*Devp + dpeq_dp*dDevp_dq - d2peq_dq2*Devpv;

    fac = (lambda-kappa)/tau*dt*(peq/peqp)^((lambda-kappa)/mu-1)/peqp;
    dpsi2_dp = dDevpv_dp - fac*dpeq_dp + fac*peq/peqp*dpeqp_dp;
    dpsi2_dq = -fac*dpeq_dq;

    J = [dpsi1_dp,dpsi1_dq;
        dpsi2_dp,dpsi2_dq];

    Dc = diag(1./sqrt(diag(J'*J)));
    J = J*Dc;

    r = [psi1;psi2];

    if (cond(J) > 1/eps)
        condIssue = true;
    else
        dx = -J\r;
        dx = Dc*dx;

        p = p + dx(1);
        q = q + dx(2);
    end

    rnorm = norm(r);
    if (iter == 1)
        rnorm0 = rnorm;
    end
    res(iter) = rnorm;

    conv = ((rnorm < rtol*rnorm0+atol) || iter >= itmax) || condIssue;

    %fprintf('%5i %15.6e\n', iter, rnorm/rnorm0);

end
flag = (rnorm < rtol*rnorm0+atol);

res = res(1:iter);
%if (iter == itmax)
%    fprintf('Non convergence with residual %15.6e\n', rnorm);
%end

% n tensor in vector form
if (qtrial > 0)
    strial = p*sbar + s0;
    nt = 3/2*strial/qtrial;
else
    nt = m;
end

% New stress
sigma = 2/3*q*nt + p*m;

phat = p + cohes/tan(fric_ang);
peq = phat + q^2/(M^2*phat);
dpeq_dp = 1 - q^2/(M^2*phat^2);
dpeq_dq = 2*q/(M^2*phat);
d2peq_dp2 = 2*q^2/(M^2*phat^3);
d2peq_dpdq = -2*q/(M^2*phat^2);
d2peq_dq2 = 2/(M^2*phat);

dptrial_dp = pbar;
dqtrial_dp = 3/(2*qtrial)*(p*sbar_sbar+sbar_s0);

B = dt*(lambda-kappa)/tau*(peq/peqp)^((lambda-kappa)/mu-1)*1/peqp;
A = 1+B*(peq/peqp)*peqp0/(lambda-kappa)*exp(Devpv/(lambda-kappa));
dlambda = Devpv / dpeq_dp;
F11 = 1 + K*dpeq_dp*B/A + Kbar*dpeq_dp*dlambda - dptrial_dp;
F12 = K*dpeq_dq*B/A;
F21 = 3*G*d2peq_dpdq*dlambda+3*G*dpeq_dq*B/A-3*G*dpeq_dq/dpeq_dp*d2peq_dp2*dlambda + 3*Gbar*dpeq_dq*dlambda - dqtrial_dp;
F22 = 1+3*G*d2peq_dq2*dlambda+3*G*dpeq_dq^2/dpeq_dp*B/A-3*G*dpeq_dq/dpeq_dp*d2peq_dpdq*dlambda;
detF = F11*F22-F12*F21;
c11 = F22/detF;
c12 = -F12/detF;
c21 = -F21/detF;
c22 = F11/detF;

edev = nt'*devp;
% sbar constains Dev*Ds*devp since sbar = 2*Gbar*Dev*Ds*devp
t = sbar/(2*Gbar);

d1 = 2/3*(c21-2*Gbar*q/qtrial*c11*edev)*K;
d2 = 2*G*c12;
d3 = 4/3*((c22-q/qtrial)-2*Gbar*q/qtrial*c12*edev)*G;
d4 = c11*K;
d5 = 2*G*q/qtrial;
d6 = 2*Gbar*q/qtrial*c11*K;
d7 = 4*Gbar*q/qtrial*c12*G;

stif = d1*(nt*m') + d2*(m*nt') + d3*(nt*nt') + d4*(m*m') + d5*Dev*Ds + d6*(t*m') + d7*(t*nt');

%stif = 0.5*(stif + stif');

end
