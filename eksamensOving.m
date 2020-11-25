A = [ 0 1; 0 1];
B = [0; 1];
C = [1 0];
o = obsv(A, C); % observabillity
r_o = rank(o); 

c = ctrb(A,B); % controlabillity
r_c= rank(c);

% K = [k_1 k_2];
% pol = eig(A-B*K); 

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
