clear, close, clc

% === Init System ===
% initial condition
x0 = [
    -1;       % position
    8*pi/7;  % angle
    0;       % velocity
    0        % angular velocity 
];

sys = CartPendulum(x0);

% === Init Control ===
% set linearization
xbar = [0, pi, 0, 0]';
ubar = 0;
[A, B] = sys.linearize_discrete(xbar, ubar);
C = eye(sys.n); % all states are measured
Q = diag([10, 1, 1, 1]);
R = 0.1;
[K, P, ~] = dlqr(A, B, Q, R);
% set control parameter
N    = 20;    % discrete horizon (Th = N*sys.samplerate)
% upper and lower bound on input and states
umin = -20 - ubar;                         % in local coordinates
umax =  20 - ubar;                         % in local coordinates
xmin = [-2, deg2rad(150), -3, -5]' - xbar; % in local coordinates
xmax = [ 2, deg2rad(210),  3,  5]' - xbar; % in local coordinates

lmin = xmin; lmax = xmax; M = eye(sys.n);  % bring in the right form for qpOASES condensing

% build parameterized condensed QP
[H, gw, gy, gu, lb, ub, lbF0, lbFw, ubF0, ubFw, F, A_bar, B_bar] = tracking_LMPC_condensing(A, B, C, Q, R, P, N, umin, umax, lmin, lmax, M);

% set refference point (global coordinates)
xref = [1, pi, 0, 0]';
uref = 0;

% === Init Simulation ===
sim_T = 3;                                    % Simulation horizon
time  = 0:sys.samplerate:sim_T;  % time steps

sys.visualize()
gif('doc/MPC.gif','DelayTime',sys.samplerate,'LoopCount',inf,'frame',gcf);

% init empty QP (for qpOASES warmstart)
QP = [];

% set constant reference signals
Xref = repmat(xref - xbar,(N+1), 1);
Uref = repmat(uref - ubar, N, 1);

% init loging arrays
U = [];
X = [];
% loop through time
for t = time
    % update parmaeterized QP
    w   = sys.measure() - xbar;
    g   = gw*w + gy*Xref + gu*Uref;
    lbF = lbF0 + lbFw*w;
    ubF = ubF0 + ubFw*w;
    % solve QP
    if isempty(QP)
        [QP,u_star,fval,exitflag,iter,lambda,auxOutput] = qpOASES_sequence( 'i',H,g,F,lb,ub,lbF,ubF );
    else
        [u_star,fval,exitflag,iter,lambda,auxOutput] = qpOASES_sequence( 'h',QP,g,lb,ub,lbF,ubF );
    end
    % unroll state evolution
    x_star = A_bar*w + B_bar*u_star;
    x_star = reshape(x_star, sys.n, N+1) + xbar;
    % get system input
    u = u_star(1) + ubar;
    % log input and measurement
    U = [U, u];
    X = [X, sys.measure()];
    % integrate and visualize
    sys.integrate(u);
    sys.visualize(x_star);
    gif;
end