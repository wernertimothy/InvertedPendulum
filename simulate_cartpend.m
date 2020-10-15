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
[A, B] = sys.linearize(xbar, ubar);
% calculate controller gain via LQR
Q = diag([10, 1, 1, 1]);
R = 0.1;
K = lqr(A, B, Q, R);
% set refference point (global coordinates)
xref = [1, pi, 0, 0]';

% === Init Simulation ===
sim_T = 3;                                    % Simulation horizon
time  = sys.samplerate:sys.samplerate:sim_T;  % time steps

sys.visualize()

% init loging arrays
U = [];
X = [];
% loop through time
for t = time
    % calculate control input
    u = -K*(sys.measure() - xref);
    % log input and measurement
    U = [U, u];
    X = [X, sys.measure()];
    % integrate and visualize
    sys.integrate(u);
    sys.visualize();
end