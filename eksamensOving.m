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
[v, lambda] = eig(A);
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
    P_nd = lyap(A,N);  
catch ME
    disp("Lyapunov not possible to evaluate")
end

%% Lyapunov - Positive definite check
try chol(P_nd)
    disp("Lyapunov positive definite")
catch ME
    disp("Lyapunov not positive definite")
end

%% LQR
Q_lqr = eye(2);
R_lqr = eye(1);
try 
    K_lqr = lqr(A, B, Q_lqr, R_lqr);
catch ME
    disp("Not possible to stabilize system using LQR");
    K_lqr = 0;
end

%% Continous KF
syms p_11 p_12 p_22;

% Variance from the given system
Q_c = 2; 
R_c = 3;

P_nd = [p_11 p_12; 
     p_12 p_22];

% Hardcoded solution of the ricatti-equation. Could be improved by using
% an in-built function
P_dot = A*P_nd + P_nd*A' + G*Q_c*G' - P_nd*C'*inv(R_c)*C*P_nd == zeros(2);
P_sol = solve(P_dot,[p_11 p_12 p_22]);
p_11 = subs(P_sol.p_11);
p_12 = subs(P_sol.p_12);
p_22 = subs(P_sol.p_22);


P_nd = [p_11 p_12;
     p_12 p_22];

% Kalman gain in the continous time
K_kf_cont = P_nd*C'*inv(R_c);


%% Discrete KF, one dimension

ntimes = 5; % Number of iterations + 1

% System matrices - given that the system is known in continous time
% Must be one-dimensional
A = -1;
B = 1; 
C = 1/10;
D = 0;

sysc = ss(A,B,C,D); % Continous state space model

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
q_d = 4; 
r_d = 25;

% Initial values
K_kf_disc_1d = zeros(1,ntimes); 
x_pri_1d = 0; 
P_pri_1d = 12; 

for k = 1:ntimes
    K_kf_disc_1d(k) = (P_pri_1d(k)*C_d')*(C_d*P_pri_1d(k)*C_d + r_d)^(-1);
    x_hat_1d(k) = x_pri_1d(k) + K_kf_disc_1d(k)*(y(k) - C_d*x_pri_1d(k));
    P_hat(k) = (eye(1) - K_kf_disc_1d(k)*C_d)*P_pri_1d(k)*(eye(1) - K_kf_disc_1d(k)*C_d)' + K_kf_disc_1d(k)*r_d*K_kf_disc_1d(k)';
    x_pri_1d(k+1) = A_d*x_hat_1d(k) + B_d*u(k);
    P_pri_1d(k+1) = A_d*P_hat(k)*A_d' + q_d;
end


%% Discrete KF, 2 or higher dimensions - DO NOT USE! DOES NOT WORK!!
%{
ntimes = 6; % Number of iterations + 1
dim = 2;    % Dimension-size
T_s = 1;    % Sample time

% System matrices 
A = [1 0; 
     0  1];
B = [0; 
     0];
C = [1 0;
     0 1];
D = [0;
     0];

% Variances in the system
q_d = 0; 
r_d = 25;

% Measurements
y_t = [10 6 1 1 11 -3;
       4 0 1 16 8 8];

% Initial values
x0 = zeros(dim,1);
x_pri_nd = x0; % Initial a priori estimate
P_pri_nd = r_d^2*eye(dim); % Initial a priori P

% Matrix for state values:
x_all = zeros(dim,ntimes);

% Discretizing 
sysc = ss(A,B,C,D);     % Continous state space model
sysd = c2d(sysc, T_s);  % Discrete system model
A_d = sysd.A;
B_d = sysd.B;
C_d = sysd.C;
D_d = sysd.D;
Q_d = q_d*eye(dim);
R_d = r_d*eye(dim);

for k = 1:ntimes
    L = P_pri_nd*C_d'*inv(C_d*P_pri_nd*C_d' + R_d);
    x_hat_nd = x_pri_nd + L*(y_t(:,k) - C_d*x_pri_nd);
    P_nd = (eye(dim) - L*C_d)*P_pri_nd*(eye(dim) - L*C_d)' + L*R_d*L';
    x_pri_nd = A_d*x_hat_nd;
    P_pri_nd = A_d*P_nd*A_d' + Q_d;
    x_all(:,k) = x_hat_nd;
end
%}


%% Continuous KF simulation

% Constants
n_states = 4;
dt = 0.01;
t_sim = 10;

t = dt:dt:t_sim;

% Continuous state space

A = zeros(n_states);
B = ones(n_states, 1);
C = eye(4);
D = zeros(size(C,1), size(B,2));

% Create system with disturbance and noise
Qw = 0.01*eye(n_states);                    % Disturbance/process noise covariance matrix
Rv = 0.1*eye(rank(C));                      % Measurement noise covariance matrix

B_dist = [B Qw 0*B];                        % System gain matrix in the noisy system
D_dist = zeros(size(C,1), size(B_dist,2));  % Feedthrough term in the noisy system
D_dist(:, end) = diag(Rv);                  % Insert noise covariance in the feedthrough matrix

ss_cont = ss(A, B_dist, C, D_dist);         % The state space with noise terms added

% Kalman filter (continuous system)
[L, P, E] = lqe(A, Qw, C, Qw, Rv);    
B_kf = [B L];

ss_kf = ss(A-L*C, B_kf, eye(n_states), 0*B_kf);

% Create simulation input
u_dist = randn(n_states, size(t, 2));
u_noise = randn(size(t));

u = 0*t; u(10:50) = 50; u(200:900) = -30;

% Augmented system input
u_aug = [u; Qw*Qw*u_dist; u_noise];

% Simulated system
[y, t] = lsim(ss_cont, u_aug, t);
plot(t, y);

