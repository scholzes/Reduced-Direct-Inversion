% clc;

% Parameters:
% n = Dimension
% N = Length of the Frame
% L = The Erasure Set

n = 250;
N = 1000;
L = [1:100];

% % The columns of F are a Gaussian randomly generated frame.
% % The columns of G are the standard dual to F.
% 
% F = (1/sqrt(n)) * randn(n,N);
% S = F * F';
% G = S \ F;

F = randn(N,n);
[F,~] = qr(F,0);
F = sqrt(N/n) * F';
G = n/N * F;

% f is a random vector that we will try to recover 
% from frame coefficient erasures.

f = rand(n,1);
f = f ./ norm(f);

% FC are the frame coefficients of f.

FC = G' * f;

% We erase the frame coefficients indexed by
% L, 

FC(L) = zeros(size(L'));

% We compute f_R.

f_R = F * FC;

% We compute the matrix M and the coefficient
% matrix.

M = (F(:,L)' * G(:,L))';
C = (eye(length(L)) - M) \ eye(length(L));

% We compute the reconstruction.

g = f_R + F(:,L) * (C * (G(:,L)' * f_R));

% We compute the \ell^2 norm of the error in the
% reconstruction.

norm(f-g)
