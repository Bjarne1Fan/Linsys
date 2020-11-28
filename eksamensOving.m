%% Matrices
A = [-1 1; 
     1  -1];
 
B = [1; 
     0];

C = [1/2 1/2];

D = 0;

% Disturbances
G = [0;
     0];
 
% Noise
H = 0;

%% Observabillity
O_obs = obsv(A, C); 
rank_o = rank(O_obs); 

%% Controlabillity
C_con = ctrb(A,B); 
rank_c = rank(C_con);

%% Eigenvalues and vectors
% Important: Matlab normalizes all eigenvectors by default. Divide by the
% lowest eigenvector to acchieve unity
[v, lambda] = eig(A)
dim = size(lambda);
%lambda = lambda * ones([dim(1),1]); % Convert from a diagonal matrix to a vector

% Finding the inverse eigenmatrix - if it exists
if det(A) ~= 0
    v_inv = inv(v)
end

%% Transfer function
sysc = ss(A,B,C,D); % Continous state space model
g = tf(sysc);

%% Inverse laplace
syms s;
time_system = ilaplace((s^2 +2)/(s^2 + 4));

%% Lyaponov
% Warning: Matlab requires unique solution
N = eye(2);
try 
    P = lyap(A,N);  
catch ME
    disp("Lyapunov not possible to evaluate")
end

%% Lyapunov - Positive definite check
try chol(P)
    disp(" ");
    disp("Lyapunov positive definite")
catch ME
    disp(" ");
    disp("Lyapunov not positive definite")
end

%% LQR
Q = eye(1);
R = eye(1);
try 
    K_lqr = lqr(A, B, Q, R);
catch ME
    disp(" ");
    disp("Not possible to stabilize system using LQR");
    K_lqr = 0;
end

%% Discrete KF
n = 5; % Number of iterations + 1

% If discretized from known system:
T_s = 1;                % Sample time
sysd = c2d(sysc, T_s);  % Discrete system model
A_d = sysd.A;
B_d = sysd.B;
C_d = sysd.C;
D_d = sysd.D;

% If system already discretized:
A_d = 1; 
B_d = 1; 
C_d = 1;
D_d = 0;

% Input and measurement
u = [100 105 92 105 102]; 
y = [10 91 220 288 405];

% Variance from the given system
Q_d = 4; 
R_d = 25;

% Initial values
K_kf_disc = zeros(1,n); 
x_pri = 0; 
P_pri = 12; 

for k = 1:n
    K_kf_disc(k) = (P_pri(k)*C_d')*(C_d*P_pri(k)*C_d + R_d)^(-1);
    x_hat(k) = x_pri(k) + K_kf_disc(k)*(y(k) - C_d*x_pri(k));
    P_hat(k) = (eye(1) - K_kf_disc(k)*C_d)*P_pri(k)*(eye(1) - K_kf_disc(k)*C_d)' + K_kf_disc(k)*R_d*K_kf_disc(k)';
    x_pri(k+1) = A_d*x_hat(k) + B_d*u(k);
    P_pri(k+1) = A_d*P_hat(k)*A_d' + Q_d;
end

%% Continous KF
syms p_11 p_12 p_22;

% Variance from the given system
Q_c = 2; 
R_c = 3;

P = [p_11 p_12; 
     p_12 p_22];

% Hardcoded solution of the ricatti-equation. Could be improved by using
% an in-built function
P_dot = A*P + P*A' + G*Q_c*G' - P*C'*inv(R_c)*C*P == zeros(2);
P_sol = solve(P_dot,[p_11 p_12 p_22]);
p_11 = subs(P_sol.p_11);
p_12 = subs(P_sol.p_12);
p_22 = subs(P_sol.p_22);


P = [p_11 p_12;
     p_12 p_22];

% Kalman gain in the continous time
K_kf_cont = P*C'*inv(R);
