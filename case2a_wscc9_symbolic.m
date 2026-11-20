% WSCC 9-bus, 3-machine algebraic observability (Case 2A)
% PMUs at generator buses 1, 2, 3. ZIP loads at buses 5, 6, 8.
% Uses the polynomial DAE model in Section V and Table I.
% Requires the Symbolic Math Toolbox.

clear; clc;

nGen = 3;
nBus = 9;
loadBus = [5 6 8];

%% 1. Variables

% machine states, derivatives and rotor-angle variables
omega = sym('omega',[nGen 1],'real');
Eqp   = sym('Eqp',  [nGen 1],'real');
Edp   = sym('Edp',  [nGen 1],'real');
Efd   = sym('Efd',  [nGen 1],'real');
VR    = sym('VR',   [nGen 1],'real');
RF    = sym('RF',   [nGen 1],'real');
ddelta = sym('ddelta',[nGen 1],'real');
domega = sym('domega',[nGen 1],'real');
dEqp   = sym('dEqp',  [nGen 1],'real');
dEdp   = sym('dEdp',  [nGen 1],'real');
dEfd   = sym('dEfd',  [nGen 1],'real');
dVR    = sym('dVR',   [nGen 1],'real');
dRF    = sym('dRF',   [nGen 1],'real');
Sdelta  = sym('Sdelta', [nGen 1],'real');
Cdelta  = sym('Cdelta', [nGen 1],'real');
dSdelta = sym('dSdelta',[nGen 1],'real');
dCdelta = sym('dCdelta',[nGen 1],'real');

% bus voltages: V, sin(theta), cos(theta) at all 9 buses
% (buses 1-3 are measured, buses 4-9 are unknown)
Vb  = sym('V',     [nBus 1],'real');
Stb = sym('Stheta',[nBus 1],'real');
Ctb = sym('Ctheta',[nBus 1],'real');

% measured generator-bus quantities and their derivatives
dV   = sym('dV',  [nGen 1],'real');  ddV   = sym('ddV',  [nGen 1],'real');
dddV = sym('dddV',[nGen 1],'real');
I    = sym('I',   [nGen 1],'real');  dI    = sym('dI',   [nGen 1],'real');
ddI  = sym('ddI', [nGen 1],'real');  dddI  = sym('dddI', [nGen 1],'real');
Sphi = sym('Sphi',[nGen 1],'real');  Cphi  = sym('Cphi', [nGen 1],'real');
dtheta   = sym('dtheta',  [nGen 1],'real');
ddtheta  = sym('ddtheta', [nGen 1],'real');
dddtheta = sym('dddtheta',[nGen 1],'real');
dphi   = sym('dphi',  [nGen 1],'real');
ddphi  = sym('ddphi', [nGen 1],'real');
dddphi = sym('dddphi',[nGen 1],'real');

% inputs (known)
TM   = sym('TM',  [nGen 1],'real');
dTM  = sym('dTM', [nGen 1],'real');
Vref = sym('Vref',[nGen 1],'real');

% machine parameters (known)
Tdo = sym('Tdo',[nGen 1],'real');  Tqo = sym('Tqo',[nGen 1],'real');
Xd  = sym('Xd', [nGen 1],'real');  Xq  = sym('Xq', [nGen 1],'real');
Xdp = sym('Xdp',[nGen 1],'real');  Xqp = sym('Xqp',[nGen 1],'real');
Rs  = sym('Rs', [nGen 1],'real');  H   = sym('H',  [nGen 1],'real');
D   = sym('D',  [nGen 1],'real');
TE  = sym('TE', [nGen 1],'real');  TF  = sym('TF', [nGen 1],'real');
TA  = sym('TA', [nGen 1],'real');  KE  = sym('KE', [nGen 1],'real');
SE  = sym('SE', [nGen 1],'real');  KF  = sym('KF', [nGen 1],'real');
KA  = sym('KA', [nGen 1],'real');
syms omega_s real

% 17 unknowns per machine and 3 voltage quantities at buses 4-9
Theta = sym([]);
for i = 1:nGen
    Theta = [Theta; omega(i); Eqp(i); Edp(i); Efd(i); VR(i); RF(i); ...
             ddelta(i); domega(i); dEqp(i); dEdp(i); dEfd(i); dVR(i); dRF(i); ...
             Sdelta(i); Cdelta(i); dSdelta(i); dCdelta(i)];
end
for j = nGen+1:nBus
    Theta = [Theta; Vb(j); Stb(j); Ctb(j)];
end

%% 2. WSCC network

% connections between WSCC buses
branch = [1 4; 2 7; 3 9; 4 5; 4 6; 5 7; 6 9; 7 8; 8 9];

% Y-bus magnitude and sine/cosine of its angle (symbolic)
adj = eye(nBus);   % 1 where Y-bus entry is nonzero
for e = 1:size(branch,1)
    adj(branch(e,1),branch(e,2)) = 1;  adj(branch(e,2),branch(e,1)) = 1;
end

Ymag = sym(zeros(nBus));  cA = sym(zeros(nBus));  sA = sym(zeros(nBus));
for b = 1:nBus
    Ymag(b,b) = sym(sprintf('Y_%d_%d',b,b),'real');
    cA(b,b)   = sym(sprintf('cA_%d_%d',b,b),'real');
    sA(b,b)   = sym(sprintf('sA_%d_%d',b,b),'real');
end
for e = 1:size(branch,1)
    b = branch(e,1);  k = branch(e,2);
    Ymag(b,k) = sym(sprintf('Y_%d_%d',b,k),'real');  Ymag(k,b) = Ymag(b,k);
    cA(b,k)   = sym(sprintf('cA_%d_%d',b,k),'real'); cA(k,b)   = cA(b,k);
    sA(b,k)   = sym(sprintf('sA_%d_%d',b,k),'real'); sA(k,b)   = sA(b,k);
end

% ZIP load parameters
PL0 = sym('PL0',[nBus 1],'real');  QL0 = sym('QL0',[nBus 1],'real');
aP  = sym('aP', [nBus 1],'real');  bP  = sym('bP', [nBus 1],'real');
cP  = sym('cP', [nBus 1],'real');
aQ  = sym('aQ', [nBus 1],'real');  bQ  = sym('bQ', [nBus 1],'real');
cQ  = sym('cQ', [nBus 1],'real');
V0  = sym('V0', [nBus 1],'real');

%% 3. Equations (all equal to zero)

P = sym([]);
stator = sym([]);

for i = 1:nGen
    % d- and q-axis current and voltage
    Id = I(i)*(Sdelta(i)*Cphi(i) - Cdelta(i)*Sphi(i));
    Iq = I(i)*(Cdelta(i)*Cphi(i) + Sdelta(i)*Sphi(i));
    Vd = Vb(i)*(Sdelta(i)*Ctb(i) - Cdelta(i)*Stb(i));
    Vq = Vb(i)*(Cdelta(i)*Ctb(i) + Sdelta(i)*Stb(i));

    % generator + exciter (7 equations)
    p1 = Tdo(i)*dEqp(i) - Efd(i) + Eqp(i) + (Xd(i)-Xdp(i))*Id;
    p2 = Tqo(i)*dEdp(i) + Edp(i) - (Xq(i)-Xqp(i))*Iq;
    p3 = ddelta(i) - omega(i) + omega_s;
    p4 = (2*H(i)/omega_s)*domega(i) - TM(i) + D(i)*(omega(i)-omega_s) ...
         + Edp(i)*Id + Eqp(i)*Iq + (Xqp(i)-Xdp(i))*Id*Iq;
    p5 = TE(i)*dEfd(i) + (KE(i)+SE(i))*Efd(i) - VR(i);
    p6 = TF(i)*dRF(i) + RF(i) - (KF(i)/TF(i))*Efd(i);
    p7 = TA(i)*dVR(i) + VR(i) - KA(i)*RF(i) + (KA(i)*KF(i)/TF(i))*Efd(i) ...
         - KA(i)*(Vref(i)-Vb(i));

    % stator (2 equations)
    p8 = Edp(i) - Vd - Rs(i)*Id + Xqp(i)*Iq;
    p9 = Eqp(i) - Vq - Rs(i)*Iq - Xdp(i)*Id;

    % sin/cos of delta (3 equations)
    p10 = Sdelta(i)^2 + Cdelta(i)^2 - 1;
    p11 = dSdelta(i) - Cdelta(i)*ddelta(i);
    p12 = dCdelta(i) + Sdelta(i)*ddelta(i);

    P = [P; p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11; p12];
    stator = [stator; p8; p9];
end

% network: P and Q balance at every bus (18 equations)
for b = 1:nBus
    sumP = sym(0);  sumQ = sym(0);
    for k = 1:nBus
        if adj(b,k)
            % cos and sin of (theta_b - theta_k - alpha_bk), Eqs. (35)-(36)
            Cbk = (Ctb(b)*Ctb(k) + Stb(b)*Stb(k))*cA(b,k) ...
                + (Stb(b)*Ctb(k) - Ctb(b)*Stb(k))*sA(b,k);
            Sbk = (Stb(b)*Ctb(k) - Ctb(b)*Stb(k))*cA(b,k) ...
                - (Ctb(b)*Ctb(k) + Stb(b)*Stb(k))*sA(b,k);
            sumP = sumP + Vb(b)*Vb(k)*Ymag(b,k)*Cbk;
            sumQ = sumQ + Vb(b)*Vb(k)*Ymag(b,k)*Sbk;
        end
    end

    % ZIP load (only at load buses)
    if ismember(b, loadBus)
        PL = PL0(b)*(aP(b)*(Vb(b)/V0(b))^2 + bP(b)*(Vb(b)/V0(b)) + cP(b));
        QL = QL0(b)*(aQ(b)*(Vb(b)/V0(b))^2 + bQ(b)*(Vb(b)/V0(b)) + cQ(b));
    else
        PL = sym(0);  QL = sym(0);
    end

    if b <= nGen
        % generator bus: injected power from the measured phasors
        Pinj = Vb(b)*I(b)*(Ctb(b)*Cphi(b) + Stb(b)*Sphi(b));
        Qinj = Vb(b)*I(b)*(Stb(b)*Cphi(b) - Ctb(b)*Sphi(b));
        P = [P; Pinj + PL - sumP; Qinj + QL - sumQ];
    else
        P = [P; PL - sumP; QL - sumQ];
    end
end

% unit circle for the unmeasured bus angles (6 equations)
for j = nGen+1:nBus
    P = [P; Stb(j)^2 + Ctb(j)^2 - 1];
end

%% 4. Stator derivatives

% use the machine equations when differentiating the stator equations
state = sym([]);  stateRate = sym([]);
for i = 1:nGen
    Id = I(i)*(Sdelta(i)*Cphi(i) - Cdelta(i)*Sphi(i));
    Iq = I(i)*(Cdelta(i)*Cphi(i) + Sdelta(i)*Sphi(i));

    fomega = (omega_s/(2*H(i)))*(TM(i) - D(i)*(omega(i)-omega_s) ...
             - Edp(i)*Id - Eqp(i)*Iq - (Xqp(i)-Xdp(i))*Id*Iq);
    fEqp = (Efd(i) - Eqp(i) - (Xd(i)-Xdp(i))*Id)/Tdo(i);
    fEdp = (-Edp(i) + (Xq(i)-Xqp(i))*Iq)/Tqo(i);
    fEfd = (VR(i) - (KE(i)+SE(i))*Efd(i))/TE(i);
    fVR  = (-VR(i) + KA(i)*RF(i) - (KA(i)*KF(i)/TF(i))*Efd(i) ...
            + KA(i)*(Vref(i)-Vb(i)))/TA(i);
    fRF  = (-RF(i) + (KF(i)/TF(i))*Efd(i))/TF(i);
    fSdelta = Cdelta(i)*(omega(i)-omega_s);
    fCdelta = -Sdelta(i)*(omega(i)-omega_s);

    state = [state; omega(i); Eqp(i); Edp(i); Efd(i); VR(i); RF(i); ...
             Sdelta(i); Cdelta(i)];
    stateRate = [stateRate; fomega; fEqp; fEdp; fEfd; fVR; fRF; ...
                 fSdelta; fCdelta];
end

% measured quantities and their time derivatives (constant Vref)
known = sym([]);  knownDerivatives = sym([]);
for i = 1:nGen
    known = [known; Vb(i); dV(i); ddV(i); I(i); dI(i); ddI(i); ...
             Stb(i); Ctb(i); Sphi(i); Cphi(i); ...
             dtheta(i); ddtheta(i); dphi(i); ddphi(i); TM(i); Vref(i)];
    knownDerivatives = [knownDerivatives; dV(i); ddV(i); dddV(i); ...
             dI(i); ddI(i); dddI(i); ...
             Ctb(i)*dtheta(i); -Stb(i)*dtheta(i); ...
             Cphi(i)*dphi(i); -Sphi(i)*dphi(i); ...
             ddtheta(i); dddtheta(i); ddphi(i); dddphi(i); dTM(i); sym(0)];
end

%% 5. Rank test

assert(length(P) == 60, 'Expected 60 initial equations.');
assert(length(Theta) == 69, 'Expected 69 latent variables.');

nLatent = length(Theta);
qMax = 3;
observable = false;

fprintf('WSCC 9-bus observability: Case 2A\n');
fprintf('Latent variables: %d (3 machines x 17 + 6 buses x 3)\n', nLatent);
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
    fprintf('The system is algebraically observable .\n');
else
    fprintf('Rank %d < %d up to q = %d.\n', r, nLatent, qMax);
    fprintf('Observability not established up to q = %d.\n', qMax);
end
