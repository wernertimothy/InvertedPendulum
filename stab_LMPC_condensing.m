%{
this is a handle for the stabilizing Linear MPC scheme.

min.   1/2 sum_{k=0}^{N-1} ||x_k||_Q + ||u_k||_R + ||x_N||_P
s.t.   x^+  == Ax + Bu        \forall k\in[0,N-1]
       lmin <= Mx_k <= lmax   \forall k\in[0,N]
       umin <= u_k  <= umax   \forall k\in[0,N-1]

The handle takes the system and controller parameter and returns
the QP parameter which describe the condensed parametric QP solved by qpOASES.

min.   1/2 U'HU + U'g(w)
s.t.   lb     <= U  <= ub
       lbF(w) <= FU <= ubF(w)

where w is the inital condition, i.e. w = x(t).
The returned problem is parameteriued in w, s.t.

        g   = gw*w
        lbF = lbF0 + lbFw*w
        ubF = ubF0 + ubFw*w

This update must be performed each time step with the new initial condition
w = x(t).
%}
function [H, gw, lb, ub, lbF0, lbFw, ubF0, ubFw, F, A_bar, B_bar] = stab_LMPC_condensing(A, B, Q, R, P, N, umin, umax, lmin, lmax, M)

[n, p] = size(B);
m      = size(M,1); 

% preallocate memory
Q_bar   = zeros((N+1)*n);
R_bar   = zeros(N*p);
A_bar   = zeros((N+1)*n, n);
B_bar   = zeros((N+1)*n, N*p);
lb      = repmat(umin, N, 1);
ub      = repmat(umax, N, 1);
ll_bar  = repmat(lmin, N+1, 1);
lu_bar  = repmat(lmax, N+1, 1);
M_bar   = zeros((N+1)*m, (N+1)*n);

% stacking the penalty
for k = 1:N
    Q_bar((k-1)*n+1:k*n, (k-1)*n+1:k*n) = Q;
    R_bar((k-1)*p+1:k*p, (k-1)*p+1:k*p) = R;
end
Q_bar(N*n+1:end, N*n+1:end) = P;

% stack constraints
for k = 1:N+1
    M_bar((k-1)*m+1:k*n, (k-1)*m+1:k*n) = M;
end

% state evolution
A_bar(1:n, 1:n) = eye(n);
tmp_B_bar = B;
for k = 1:N
    A_bar(k*n+1:(k+1)*n, : )  = A^k;
    B_bar(k*n+1:(k+1)*n, 1:k*p) = tmp_B_bar;
    tmp_B_bar = [A^k*B, tmp_B_bar];
end

% constraints
lbF0 = ll_bar;
lbFw = -M_bar*A_bar;
ubF0 = lu_bar;
ubFw = -M_bar*A_bar;
F    = B_bar;

% calculate condensed problem
H   = B_bar'*Q_bar*B_bar + R_bar;
gw  = B_bar'*Q_bar*A_bar;
end
