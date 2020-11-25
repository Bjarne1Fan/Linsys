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
lambda = lambda * ones([dim(1),1]);

%% transferfunction
A = [0 -4; 1 0]; B = [-2; 0]; C = [0 1]; D = 1;
sys = ss(A,B,C,D); % state space model
g = tf(sys);

%% inverse laplace
syms s;
re = ilaplace((s^2 +2)/(s^2 + 4));

%% positive definite check
M = [1 0; 0 1];
try chol(M)
    disp("Positive definite")
catch ME
    disp("Not positve definite")
end
