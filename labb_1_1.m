
kappa = 2*pi;
N=6;
M=10000;
x_ps = unifrnd(0,1,M,2);


% the actual v in the assignment
v=@(x_p) x_p*[[0,-1];[1,0]] + [0,1];
    
% to test with the characteristic we calculated in (a)
% v(x_p) = x_p
%v = @(x_p) x_p;

f=@(x_p) double(vecnorm(x_p-0.5,2,2)<0.1);

v_ps=zeros(M,2);
for i=1:M
    v_ps(i,:) = v(x_ps(i,:));
end

f_ps=f(x_ps);

% we fill the matrix using a double loop
width = (2*N+1);
A = zeros(M, width^2);
index = @(n_1,n_2) width*(N+n_1)+(n_2+N)+1; % makes number from {-N,...,N} x {-N,...,N} to {1,...,(2N+1)^2}


for m = 1:M
    for n_1 = -N : N
        for n_2 = -N : N
            A(m, index(n_1,n_2)) = exp(kappa*1i*(n_1 * x_ps(m,1) + n_2 * x_ps(m,2)));
        end
    end
end

v_ps_1 = v_ps(:,1);
v_ps_2 = v_ps(:,2);

% the coeffs come in the wrong order for some reason
coeffs_1 = flip(A \ v_ps_1);
coeffs_2 = flip(A \ v_ps_2);
coeffs_f = flip(A \ f_ps);

exponent_vector = zeros(width*width,2);
    for n_1 = -N : N
        for n_2 = -N : N
            exponent_vector(index(n_1,n_2),:) = 1i * kappa * [n_1, n_2];
        end
    end

e = exp(1);
v_1 = @(x,y) real(coeffs_1' * e .^ (exponent_vector * [x;y]));
v_2 = @(x,y) real(coeffs_2' * e .^ (exponent_vector * [x;y]));

v_numerical = @(x,y) [v_1(x,y); v_2(x,y)];
f_numerical = @(x,y) real(coeffs_f' * e.^ (exponent_vector * [x;y]));

%% We plot v and v_numerical to make sure we have made calculated correctly
S=15; % we take some epsilon off the corner because fourier interpolation is scuffed on edges
x = linspace(0.1,1-0.1,S);
y = linspace(0.1,1-0.1,S);
[X, Y] = meshgrid(x,y);

Z1=zeros(S,S,2);
Z2=zeros(S,S,2);
for i=1:S
    for j=1:S
        Z1(i,j,:) = v([X(i,j), Y(i,j)]);
        Z2(i,j,:) = (v_numerical(X(i,j), Y(i,j)));
    end
end

figure(1) % true plot
quiver(X,Y, Z1(:,:,1), Z1(:,:,2))
figure(2) % least squares interpolated plot
quiver(X,Y, Z2(:,:,1), Z2(:,:,2))


%% Now to solve for density u
K=10;
h = 1/K;
charX = zeros(K);
ks = 1:K;
[Ks1, Ks2] = meshgrid(ks,ks);
ks = [Ks1(:), Ks2(:)];
charX 

for i=1:K
    charX(i) = 
end