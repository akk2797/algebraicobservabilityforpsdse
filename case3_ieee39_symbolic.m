% IEEE 39-bus, 10-machine algebraic observability (Case 3)
% PMUs at all generator buses (30-39).
% Network and loads are taken from MATPOWER case39. No Kron reduction.
% Requires Symbolic Math Toolbox and case39.m on the MATLAB path.

clear; clc;

%% 1. Variables

nBus = 39;
nGen = 10;

% Machine states, first derivatives, and rotor-angle variables
omega  = sym('omega',  [nGen 1], 'real');
Eqp    = sym('Eqp',    [nGen 1], 'real');
Edp    = sym('Edp',    [nGen 1], 'real');
Efd    = sym('Efd',    [nGen 1], 'real');
VR     = sym('VR',     [nGen 1], 'real');
RF     = sym('RF',     [nGen 1], 'real');
ddelta = sym('ddelta', [nGen 1], 'real');
domega = sym('domega', [nGen 1], 'real');
dEqp   = sym('dEqp',   [nGen 1], 'real');
dEdp   = sym('dEdp',   [nGen 1], 'real');
dEfd   = sym('dEfd',   [nGen 1], 'real');
dVR    = sym('dVR',    [nGen 1], 'real');
dRF    = sym('dRF',    [nGen 1], 'real');
Sdelta  = sym('Sdelta',  [nGen 1], 'real');
Cdelta  = sym('Cdelta',  [nGen 1], 'real');
dSdelta = sym('dSdelta', [nGen 1], 'real');
dCdelta = sym('dCdelta', [nGen 1], 'real');

% Bus voltages and voltage-angle sine/cosine
Vb  = sym('V',      [nBus 1], 'real');
Stb = sym('Stheta', [nBus 1], 'real');
Ctb = sym('Ctheta', [nBus 1], 'real');

% Measured current phasors and their derivatives
I     = sym('I',    [nGen 1], 'real');
dV    = sym('dV',   [nGen 1], 'real');
ddV   = sym('ddV',  [nGen 1], 'real');
dddV  = sym('dddV', [nGen 1], 'real');
dI    = sym('dI',   [nGen 1], 'real');
ddI   = sym('ddI',  [nGen 1], 'real');
dddI  = sym('dddI', [nGen 1], 'real');
Sphi  = sym('Sphi', [nGen 1], 'real');
Cphi  = sym('Cphi', [nGen 1], 'real');
dtheta   = sym('dtheta',   [nGen 1], 'real');
ddtheta  = sym('ddtheta',  [nGen 1], 'real');
dddtheta = sym('dddtheta', [nGen 1], 'real');
dphi     = sym('dphi',     [nGen 1], 'real');
ddphi    = sym('ddphi',    [nGen 1], 'real');
dddphi   = sym('dddphi',   [nGen 1], 'real');

% Known inputs and machine parameters
TM   = sym('TM',   [nGen 1], 'real');
dTM  = sym('dTM',  [nGen 1], 'real');
Vref = sym('Vref', [nGen 1], 'real');
Tdo = sym('Tdo', [nGen 1], 'real');  Tqo = sym('Tqo', [nGen 1], 'real');
Xd  = sym('Xd',  [nGen 1], 'real');  Xq  = sym('Xq',  [nGen 1], 'real');
Xdp = sym('Xdp', [nGen 1], 'real');  Xqp = sym('Xqp', [nGen 1], 'real');
Rs  = sym('Rs',  [nGen 1], 'real');  H   = sym('H',   [nGen 1], 'real');
D   = sym('D',   [nGen 1], 'real');
TE  = sym('TE',  [nGen 1], 'real');  TF  = sym('TF',  [nGen 1], 'real');
TA  = sym('TA',  [nGen 1], 'real');  KE  = sym('KE',  [nGen 1], 'real');
SE  = sym('SE',  [nGen 1], 'real');  KF  = sym('KF',  [nGen 1], 'real');
KA  = sym('KA',  [nGen 1], 'real');
syms omega_s real

%% 2. IEEE 39-bus network

% Only the case39 data file is needed; no MATPOWER functions are called.
assert(exist('case39', 'file') == 2, ...
    'Place the MATPOWER case39.m file in the current folder or MATLAB path.');
mpc = case39;
assert(size(mpc.bus,1) == nBus && size(mpc.gen,1) == nGen, ...
    'The case file must have 39 buses and 10 generators.');

baseMVA = mpc.baseMVA;
busNumber = mpc.bus(:,1);
assert(isequal(busNumber, (1:nBus).'), 'Bus numbers must be 1 through 39.');
genBus = mpc.gen(:,1);
assert(numel(unique(genBus)) == nGen, 'Expected one machine per generator bus.');
assert(isequal(sort(genBus), (30:39).'), 'Expected generator buses 30 through 39.');
otherBus = setdiff((1:nBus).', genBus);

% Build Ybus from line impedances, taps, charging, and bus shunts
Ybus = complex(zeros(nBus,nBus));
for e = 1:size(mpc.branch,1)
    if mpc.branch(e,11) == 0
        continue
    end
    b = mpc.branch(e,1);
    k = mpc.branch(e,2);
    y = 1/(mpc.branch(e,3) + 1i*mpc.branch(e,4));
    ysh = 1i*mpc.branch(e,5)/2;
    tap = mpc.branch(e,9);
    if tap == 0
        tap = 1;
    end
    tap = tap*exp(1i*pi/180*mpc.branch(e,10));
    Ybus(b,b) = Ybus(b,b) + (y+ysh)/abs(tap)^2;
    Ybus(b,k) = Ybus(b,k) - y/conj(tap);
    Ybus(k,b) = Ybus(k,b) - y/tap;
    Ybus(k,k) = Ybus(k,k) + y+ysh;
end
for b = 1:nBus
    Ybus(b,b) = Ybus(b,b) + (mpc.bus(b,5)+1i*mpc.bus(b,6))/baseMVA;
end

% Network coefficients are fixed; all unknown variables remain symbolic
G = real(Ybus);
B = imag(Ybus);

% Signed load injection in per unit (negative for consumption)
PL = -mpc.bus(:,3)/baseMVA;
QL = -mpc.bus(:,4)/baseMVA;

% 17 latent variables per generator and 3 at each unmeasured bus
Theta = sym([]);
for i = 1:nGen
    Theta = [Theta; omega(i); Eqp(i); Edp(i); Efd(i); VR(i); RF(i); ...
             ddelta(i); domega(i); dEqp(i); dEdp(i); dEfd(i); dVR(i); dRF(i); ...
             Sdelta(i); Cdelta(i); dSdelta(i); dCdelta(i)];
end
for j = otherBus.'
    Theta = [Theta; Vb(j); Stb(j); Ctb(j)];
end

%% 3. Equations (all equal to zero)

P = sym([]);
stator = sym([]);

for i = 1:nGen
    b = genBus(i);

    Id = I(i)*(Sdelta(i)*Cphi(i) - Cdelta(i)*Sphi(i));
    Iq = I(i)*(Cdelta(i)*Cphi(i) + Sdelta(i)*Sphi(i));
    Vd = Vb(b)*(Sdelta(i)*Ctb(b) - Cdelta(i)*Stb(b));
    Vq = Vb(b)*(Cdelta(i)*Ctb(b) + Sdelta(i)*Stb(b));

    % Generator and exciter equations
    p1 = Tdo(i)*dEqp(i) - Efd(i) + Eqp(i) + (Xd(i)-Xdp(i))*Id;
    p2 = Tqo(i)*dEdp(i) + Edp(i) - (Xq(i)-Xqp(i))*Iq;
    p3 = ddelta(i) - omega(i) + omega_s;
    p4 = (2*H(i)/omega_s)*domega(i) - TM(i) + D(i)*(omega(i)-omega_s) ...
         + Edp(i)*Id + Eqp(i)*Iq + (Xqp(i)-Xdp(i))*Id*Iq;
    p5 = TE(i)*dEfd(i) + (KE(i)+SE(i))*Efd(i) - VR(i);
    p6 = TF(i)*dRF(i) + RF(i) - (KF(i)/TF(i))*Efd(i);
    p7 = TA(i)*dVR(i) + VR(i) - KA(i)*RF(i) + (KA(i)*KF(i)/TF(i))*Efd(i) ...
         - KA(i)*(Vref(i)-Vb(b));

    % Stator equations
    p8 = Edp(i) - Vd - Rs(i)*Id + Xqp(i)*Iq;
    p9 = Eqp(i) - Vq - Rs(i)*Iq - Xdp(i)*Id;

    % Rotor-angle equations
    p10 = Sdelta(i)^2 + Cdelta(i)^2 - 1;
    p11 = dSdelta(i) - Cdelta(i)*ddelta(i);
    p12 = dCdelta(i) + Sdelta(i)*ddelta(i);

    P = [P; p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11; p12];
    stator = [stator; p8; p9];
end

% Active- and reactive-power balance at all 39 buses
for b = 1:nBus
    sumP = sym(0);
    sumQ = sym(0);
    neighbor = find(Ybus(b,:) ~= 0);
    for k = neighbor
        Cdiff = Ctb(b)*Ctb(k) + Stb(b)*Stb(k);
        Sdiff = Stb(b)*Ctb(k) - Ctb(b)*Stb(k);
        sumP = sumP + Vb(b)*Vb(k)*(G(b,k)*Cdiff + B(b,k)*Sdiff);
        sumQ = sumQ + Vb(b)*Vb(k)*(G(b,k)*Sdiff - B(b,k)*Cdiff);
    end

    i = find(genBus == b, 1);
    if ~isempty(i)
        Pinj = Vb(b)*I(i)*(Ctb(b)*Cphi(i) + Stb(b)*Sphi(i));
        Qinj = Vb(b)*I(i)*(Stb(b)*Cphi(i) - Ctb(b)*Sphi(i));
        P = [P; Pinj + PL(b) - sumP; Qinj + QL(b) - sumQ];
    else
        P = [P; PL(b) - sumP; QL(b) - sumQ];
    end
end

% Unit circle for voltage angles at unmeasured buses
for j = otherBus.'
    P = [P; Stb(j)^2 + Ctb(j)^2 - 1];
end

%% 4. Stator derivatives

% Use the machine equations to eliminate higher state derivatives
state = sym([]);
stateRate = sym([]);
for i = 1:nGen
    b = genBus(i);
    Id = I(i)*(Sdelta(i)*Cphi(i) - Cdelta(i)*Sphi(i));
    Iq = I(i)*(Cdelta(i)*Cphi(i) + Sdelta(i)*Sphi(i));

    fomega = (omega_s/(2*H(i)))*(TM(i) - D(i)*(omega(i)-omega_s) ...
             - Edp(i)*Id - Eqp(i)*Iq - (Xqp(i)-Xdp(i))*Id*Iq);
    fEqp = (Efd(i) - Eqp(i) - (Xd(i)-Xdp(i))*Id)/Tdo(i);
    fEdp = (-Edp(i) + (Xq(i)-Xqp(i))*Iq)/Tqo(i);
    fEfd = (VR(i) - (KE(i)+SE(i))*Efd(i))/TE(i);
    fVR = (-VR(i) + KA(i)*RF(i) - (KA(i)*KF(i)/TF(i))*Efd(i) ...
          + KA(i)*(Vref(i)-Vb(b)))/TA(i);
    fRF = (-RF(i) + (KF(i)/TF(i))*Efd(i))/TF(i);
    fSdelta = Cdelta(i)*(omega(i)-omega_s);
    fCdelta = -Sdelta(i)*(omega(i)-omega_s);

    state = [state; omega(i); Eqp(i); Edp(i); Efd(i); VR(i); RF(i); ...
             Sdelta(i); Cdelta(i)];
    stateRate = [stateRate; fomega; fEqp; fEdp; fEfd; fVR; fRF; ...
                 fSdelta; fCdelta];
end

% All generator-bus voltage/current phasors are measured
% Voltage references are treated as constant
known = sym([]);
knownDerivatives = sym([]);
for i = 1:nGen
    b = genBus(i);
    known = [known; Vb(b); dV(i); ddV(i); I(i); dI(i); ddI(i); ...
             Stb(b); Ctb(b); Sphi(i); Cphi(i); ...
             dtheta(i); ddtheta(i); dphi(i); ddphi(i); TM(i); Vref(i)];
    knownDerivatives = [knownDerivatives; dV(i); ddV(i); dddV(i); ...
             dI(i); ddI(i); dddI(i); ...
             Ctb(b)*dtheta(i); -Stb(b)*dtheta(i); ...
             Cphi(i)*dphi(i); -Sphi(i)*dphi(i); ...
             ddtheta(i); dddtheta(i); ddphi(i); dddphi(i); dTM(i); sym(0)];
end

%% 5. Symbolic rank test

assert(length(P) == 227, 'Expected 227 initial equations.');
assert(length(Theta) == 257, 'Expected 257 latent variables.');

nLatent = length(Theta);
qMax = 3;
observable = false;

fprintf('IEEE 39-bus observability: Case 3\n');
fprintf('PMUs at all 10 generator buses\n');
fprintf('Latent variables: %d (10 machines x 17 + 29 buses x 3)\n', nLatent);
fprintf('Rank needed for observability: %d\n\n', nLatent);

for q = 0:qMax
    fprintf('Computing symbolic rank at q = %d ...\n', q);
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

    % Differentiate only the 20 stator equations
    stator = jacobian(stator, state)*stateRate ...
           + jacobian(stator, known)*knownDerivatives;
    P = [P; stator];
end

fprintf('\n');
if observable
    fprintf('Full rank (%d) reached at q = %d.\n', nLatent, q);
    fprintf('The system is algebraically observable (generic full rank).\n');
else
    fprintf('Rank %d < %d up to q = %d.\n', r, nLatent, qMax);
    fprintf('Observability not established up to q = %d.\n', qMax);
end
