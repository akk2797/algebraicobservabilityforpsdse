% Single-machine algebraic observability (Case 1)
% Paper: "Algebraic Observability for Power System Dynamic State Estimation"
% Model from Section V and Table I.
%
% 4th-order generator + IEEE Type-1 exciter, two stator equations.
% Known: V, theta, I, phi (PMU at the terminal), T_M and V_ref.
% sin(delta), cos(delta) are used in place of delta so that
% everything is polynomial.
% Only the stator equations are differentiated. Machine derivatives
% that show up are replaced using the machine equations, so the
% number of unknowns stays at 17.
%
% Needs the Symbolic Math Toolbox.

clear; clc;

%% 1. Variables

% unknowns: machine states, their first derivatives, rotor-angle lift
syms omega Eqp Edp Efd VR RF real
syms ddelta domega dEqp dEdp dEfd dVR dRF real
syms Sdelta Cdelta dSdelta dCdelta real

% PMU measurements and their derivatives (known)
syms V dV ddV dddV I dI ddI dddI real
syms Stheta Ctheta Sphi Cphi real
syms dtheta ddtheta dddtheta dphi ddphi dddphi real

% inputs and machine parameters (known)
syms TM dTM Vref real
syms Tdo Tqo Xd Xq Xdp Xqp Rs H D omega_s real
syms TE TF TA KE SE KF KA real

% the 17 unknowns, same order as Eq. (37)
Theta = [omega; Eqp; Edp; Efd; VR; RF; ...
         ddelta; domega; dEqp; dEdp; dEfd; dVR; dRF; ...
         Sdelta; Cdelta; dSdelta; dCdelta];

%% 2. Equations (all written as ... = 0)

% d- and q-axis current and voltage
Id = I*(Sdelta*Cphi - Cdelta*Sphi);
Iq = I*(Cdelta*Cphi + Sdelta*Sphi);
Vd = V*(Sdelta*Ctheta - Cdelta*Stheta);
Vq = V*(Cdelta*Ctheta + Sdelta*Stheta);

% generator + exciter (7 equations)
p1 = Tdo*dEqp - Efd + Eqp + (Xd-Xdp)*Id;
p2 = Tqo*dEdp + Edp - (Xq-Xqp)*Iq;
p3 = ddelta - omega + omega_s;
p4 = (2*H/omega_s)*domega - TM + D*(omega-omega_s) ...
     + Edp*Id + Eqp*Iq + (Xqp-Xdp)*Id*Iq;
p5 = TE*dEfd + (KE+SE)*Efd - VR;
p6 = TF*dRF + RF - (KF/TF)*Efd;
p7 = TA*dVR + VR - KA*RF + (KA*KF/TF)*Efd - KA*(Vref-V);

% stator (2 equations)
p8 = Edp - Vd - Rs*Id + Xqp*Iq;
p9 = Eqp - Vq - Rs*Iq - Xdp*Id;

% sin/cos of delta: unit circle and their derivatives (3 equations)
p10 = Sdelta^2 + Cdelta^2 - 1;
p11 = dSdelta - Cdelta*ddelta;
p12 = dCdelta + Sdelta*ddelta;

% starting set P[0], 12 equations
P = [p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11; p12];
stator = [p8; p9];

%% 3. Differentiating the stator equations

% right-hand sides of the machine equations, used to replace
% state derivatives when we differentiate the stator equations
fomega = (omega_s/(2*H))*(TM - D*(omega-omega_s) ...
         - Edp*Id - Eqp*Iq - (Xqp-Xdp)*Id*Iq);
fEqp = (Efd - Eqp - (Xd-Xdp)*Id)/Tdo;
fEdp = (-Edp + (Xq-Xqp)*Iq)/Tqo;
fEfd = (VR - (KE+SE)*Efd)/TE;
fVR  = (-VR + KA*RF - (KA*KF/TF)*Efd + KA*(Vref-V))/TA;
fRF  = (-RF + (KF/TF)*Efd)/TF;
fSdelta = Cdelta*(omega-omega_s);
fCdelta = -Sdelta*(omega-omega_s);

state     = [omega; Eqp; Edp; Efd; VR; RF; Sdelta; Cdelta];
stateRate = [fomega; fEqp; fEdp; fEfd; fVR; fRF; fSdelta; fCdelta];

% known quantities and what their time derivative is
% (V_ref is taken constant)
known = [V; dV; ddV; I; dI; ddI; ...
         Stheta; Ctheta; Sphi; Cphi; ...
         dtheta; ddtheta; dphi; ddphi; TM; Vref];
knownDerivatives = [dV; ddV; dddV; dI; ddI; dddI; ...
                    Ctheta*dtheta; -Stheta*dtheta; ...
                    Cphi*dphi; -Sphi*dphi; ...
                    ddtheta; dddtheta; ddphi; dddphi; dTM; sym(0)];

%% 4. Rank test

nLatent = length(Theta);   % rank we need to reach
qMax = 3;
observable = false;

fprintf('Latent variables: %d\n', nLatent);
fprintf('Rank needed for observability: %d\n\n', nLatent);

for q = 0:qMax
    J = jacobian(P, Theta);
    r = rank(J);
    fprintf('q = %d: %d equations, %d unknowns, rank = %d, missing = %d\n', ...
            q, length(P), nLatent, r, nLatent - r);

    if r == nLatent
        observable = true;
        break
    end
    if q == qMax
        break
    end

    % take one more time derivative of the stator equations
    stator = jacobian(stator, state)*stateRate ...
           + jacobian(stator, known)*knownDerivatives;
    P = [P; stator];
end

fprintf('\n');
if observable
    fprintf('Full rank (%d) reached at q = %d.\n', nLatent, q);
    fprintf('The system is algebraically observable.\n');
else
    fprintf('Rank %d < %d up to q = %d.\n', r, nLatent, qMax);
    fprintf('The system is not algebraically observable for q <= %d.\n', qMax);
end