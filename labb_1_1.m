
kappa = 2*pi;
N=6;
M=1000;

x_ps = unifrnd(0,1,2,M);

% the actual v in the assignment
%v=@(x_p) [[0,1];[-1,0]]*x_p + [0;1];
    
% to test with the characteristic we calculated in (a)
% v(x_p) = x_p
v = @(x_p) x_p;

f=@(x_p) double(vecnorm(x_p-0.5,2,1)<0.1);

v_ps=v(x_ps);
f_ps=f(x_ps);

% we fill the matrix using a double loop
width = (2*N+1);
A = zeros(M, width^2);
index = @(n_1,n_2) width*(N+n_1)+(n_2+N)+1; % makes number from {-N,...,N} x {-N,...,N} to {1,...,(2N+1)^2}


for m = 1:M
    for n_1 = -N : N
        for n_2 = -N : N
            A(m, index(n_1,n_2)) = exp(kappa*1i*(n_1 * x_ps(1,m) + n_2 * x_ps(2,m)));
        end
    end
end

v_ps_1 = v_ps(1,:)';
v_ps_2 = v_ps(2,:)';

coeffs_1 = A \ v_ps_1;
coeffs_2 = A \ v_ps_2;
coeffs_f = A \ f_ps';

tmp = zeros(width^2,2); % black magic to create v functions
for n_1 = -N : N
    for n_2 = -N : N
        tmp(index(n_1,n_2), :) = [n_1,n_2];
    end
end
e = exp(1);
v_1 = @(x,y) real(coeffs_1' * e .^ (1i*kappa* tmp * [x;y]));
v_2 = @(x,y) real(coeffs_2' * e .^ (1i*kappa* tmp * [x;y]));

v_numerical = @(x,y) [v_1(x,y); v_2(x,y)];


v([0.5;0.70])
v_numerical(0.5,0.70)


