clear, close, clc

% Init Cartpend
x0 = [
    0;
    7*pi/8;
    0;
    0
];

sys = CartPendulum(x0);
% [A, B] = sys.linearize([0, pi-0.01, 0, 0.1], 0);
% Q = diag([1, 10, 1, 100]);
% R = 0.1;
% K = lqr(A, B, Q, R);


sim_T = 20;
time  = sys.samplerate:sys.samplerate:sim_T;

X = x0;
sys.visualize()

for t = time
%     u = -K*(sys.measure()-[0, pi, 0, 0]');
    u = 0;
    sys.integrate(u);
    sys.visualize();
end