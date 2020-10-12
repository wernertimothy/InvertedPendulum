clear, close, clc

% Init Cartpend
x0 = [
    1;
    7*pi/8;
    0;
    0
];

sys = CartPendulum(x0);
xbar = [-1, pi, 0, 0]';
ubar = 0;
[A, B] = sys.linearize(xbar, ubar);
Q = diag([10, 1, 1, 1]);
R = 0.1;
K = lqr(A, B, Q, R);


sim_T = 3;
time  = sys.samplerate:sys.samplerate:sim_T;

X = x0;
sys.visualize()

for t = time
    u = -K*(sys.measure()-xbar);
%     u = 0;
    sys.integrate(u);
    sys.visualize();
end