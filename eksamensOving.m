%% Matrices
A = [1 0; 
     0 1];
 
B = [0; 
     1];

C = [1 0];

D = 0;

% Disturbances
G = [0;
     0];
 
% Noise
H = 0;

%% Observabillity
o = obsv(A, C); 
r_o = rank(o); 

%% Controlabillity
c = ctrb(A,B); 
r_c = rank(c);

%% Eigenvalues and vectors
% Important: Matlab normalizes all eigenvectors by default. Divide by the
% lowest eigenvector to acchieve unity
[v, lambda] = eig(A);
dim = size(lambda);
lambda = lambda * ones([dim(1),1]); % Convert from a diagonal matrix to a vector

%% Transfer function
sysc = ss(A,B,C,D); % Continous state space model
g = tf(sysc);

%% Inverse laplace
syms s;
re = ilaplace((s^2 +2)/(s^2 + 4));

%% Lyaponov
% Warning: Matlab requires unique solution
N = eye(2);
M = lyap(A,N);  

%% Positive definite check
try chol(M)
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
n = 3; %number of iterations - 1

% If discretized from known system:
T_s = 1;    % Sample time
sysd = c2d(sysc, T_s);
A_d = sysd.A;
B_d = sysd.B;
C_d = sysd.C;
D_d = sysd.D;

% If not known system:
A_d = 1; 
B_d = 1; 
C_d = 1;
D_d = 0;

% Input and measurement
u = zeros(1,n); 
y = [1 2 1];

% Variance from the given system
Q_d = 2; 
R_d = 3;

% Initial values
K_kf = zeros(1,n); 
x_bar = 0; 
P_bar = eye(1); 

for k = 1:n
    K_kf(k) = (P_bar(k)*C_d')*(C_d*P_bar(k)*C_d + R_d)^(-1);
    x_hat(k) = x_bar(k) + K_kf(k)*(y(k) - C_d*x_bar(k));
    P_hat(k) = (eye(1) -K_kf(k)*C_d)*P_bar(k)*(eye(1)-K_kf(k)*C_d)' + K_kf(k)*R_d*K_kf(k)';
    x_bar(k+1) = A_d*x_hat(k) + B_d*u(k);
    P_bar(k+1) = A_d*P_hat(k)*A_d' + Q_d;
end

%% Other
syms m_1 m_2 m_3;
 
M = [m_1 m_2; 
     m_2 m_3];
M_merket = M^(-1);
