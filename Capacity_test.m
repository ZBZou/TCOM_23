clear all;
close all;

x1 = [1 -1];
x2 = [0.707+0.707i, 0.707-0.707i, -0.707+0.707i, -0.707-0.707i];
M = [1, 1i; -1i, 1];
% X = [ones(1,4),-1*ones(1,4);x2,x2];
X = [1*ones(1,4),-1*ones(1,4),1i*ones(1,4),-1i*ones(1,4);x2,x2,x2,x2];
r = M*X;
Q = reshape(real(r),[],1);
I = reshape(imag(r),[],1);
g = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16];
gscatter(Q,I,g,'rkgbcm','+oxs<');

% r1 = [1, -1, 1i, -1i];
% r2 = [-1i, 1i, 1, -1];
% x = [1, 0; -1,0; 0,1; 0,-1]

% h = 1+1i;
% [U,Sig,V] = svd([1, -1; 1, 1]);
% sig = diag(Sig);
% Phi = V(1,:) +1i*V(2,:);
% Psi = U(1,:) +1i*U(2,:);
% tx = x1(1)*Phi(1)/sig(1)+x2(1)*Phi(2)/sig(2)
% rx = h*tx;
% hat_x = Psi'*rx;