%{
this is a handle for the stabilizing Linear MPC scheme.

min.   1/2 sum_{k=0}^{N-1} ||y_k - y_r||_Q + ||u_k - u_r||_R + ||x_N - x_r||_P
s.t.   x^+  == Ax + Bu        \forall k\in[0,N-1]
       lmin <= Mx_k <= lmax   \forall k\in[0,N]
       umin <= u_k  <= umax   \forall k\in[0,N-1]

The handle takes the system and controller parameter and returns
the QP parameter which describe the condensed parametric QP solved by qpOASES.

min.   1/2 U'HU + U'g(w, yr, ur)
s.t.   lb     <= U  <= ub
       lbF(w) <= FU <= ubF(w)

where w is the inital condition, i.e. w = x(t).
The returned problem is parameteriued in w, s.t.

        g   = gw*w + gy*yr + gu*ur
        lbF = lbF0 + lbFw*w
        ubF = ubF0 + ubFw*w

This update must be performed each time step with the new initial condition
w = x(t).
%}
function [H, gw, gy, gu, lb, ub, lbF0, lbFw, ubF0, ubFw, F, A_bar, B_bar] = tracking_LMPC_condensing(A, B, C, Q, R, P, N, umin, umax, lmin, lmax, M)

[n, p] = size(B);
m      = size(M,1);
q      = size(C,1);

% preallocate memory
Q_bar   = zeros((N+1)*n);
R_bar   = zeros(N*p);
A_bar   = zeros((N+1)*n, n);
B_bar   = zeros((N+1)*n, N*p);
T_bar   = zeros((N+1)*n, (N+1)*q);
lb      = repmat(umin, N, 1);
ub      = repmat(umax, N, 1);
ll_bar  = repmat(lmin, N+1, 1);
lu_bar  = repmat(lmax, N+1, 1);
M_bar   = zeros((N+1)*m, (N+1)*n);

% stack matricies by looping through the stages 0:N-1
A_bar(1:n, 1:n) = eye(n);
tmp_B_bar = B;
for stage = 1:N
    % stack penalty
    Q_bar((stage-1)*n+1:stage*n, (stage-1)*n+1:stage*n) = C'*Q*C;
    R_bar((stage-1)*p+1:stage*p, (stage-1)*p+1:stage*p) = R;
    T_bar((stage-1)*n+1:stage*n, (stage-1)*q+1:stage*q) = C'*Q;
    % stack constraints
    M_bar((stage-1)*m+1:stage*n, (stage-1)*m+1:stage*n) = M;
    % state evolution
    A_bar(stage*n+1:(stage+1)*n, : )  = A^stage;
    B_bar(stage*n+1:(stage+1)*n, 1:stage*p) = tmp_B_bar;
    tmp_B_bar = [A^stage*B, tmp_B_bar];
end
% stage N
Q_bar(N*n+1:end, N*n+1:end) = C'*P*C;
T_bar(N*n+1:end, N*q+1:end) = C'*P;
M_bar(N*m+1:end, N*m+1:end) = M;

% constraints
lbF0 = ll_bar;
lbFw = -M_bar*A_bar;
ubF0 = lu_bar;
ubFw = -M_bar*A_bar;
F    = B_bar;

% calculate condensed problem
H   = B_bar'*Q_bar*B_bar + R_bar;
gw  = B_bar'*Q_bar*A_bar;
gy  = -B_bar'*T_bar;
gu  = -R_bar;
end
