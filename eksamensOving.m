A = [ 0 1; 0 1];
B = [0; 1];
C = [1 0];

%% observabillity
o = obsv(A, C); 
r_o = rank(o); 

%% controlabillity
c = ctrb(A,B); 
r_c= rank(c);

%% Eigenvalues and vectors
[v, lambda ]= eig(A);
dim = size(lambda);
lambda = lambda * ones([dim(1),1]); % convert from a diagnoal matrix to a vector

%% transferfunction
A = [0 -4; 1 0]; B = [-2; 0]; C = [0 1]; D = 1;
sys = ss(A,B,C,D); % state space model
g = tf(sys);

%% inverse laplace
syms s;
re = ilaplace((s^2 +2)/(s^2 + 4));

%% lyaponov
% N = eye(2);
% M = lyap(A,N);  %requires unike solution

%% positive definite check
% M = [1 0; 0 1];
% try chol(M)
%     disp("Positive definite")
% catch ME
%     disp("Not positve definite")
% end
%% LQR

%% discrete KF
n = 3; %number of iterations - 1
A_d = 1; B_d = 1; C_d = 1;
u = zeros(1,n); y = [1 2 1];
Q_d = 2; R_d = 3;
K=zeros(1,n); x_bar = 0; P_bar = eye(1); %inital values
for k = 1: n
    K(k) = (P_bar(k)*C_d')*(C_d*P_bar(k)*C_d + R_d)^(-1);
    x_hat(k) = x_bar(k) + K(k)*(y(k) - C_d*x_bar(k));
    P_hat(k) = (eye(1) -K(k)*C_d)*P_bar(k)*(eye(1)-K(k)*C_d)' + K(k)*R_d*K(k)';
    x_bar(k+1) = A_d*x_hat(k) + B_d*u(k);
    P_bar(k+1) = A_d*P_hat(k)*A_d' + Q_d;
end

%% Minimal realization



