function [X,L] = DPC(s,H)
N = size(H,1);
[Q,R] = qr(H');
L=R'; Q=Q';

xp = s;
for i = 2:N
    xp(i) = xp(i) - L(i,1:i-1)/L(i,i)*xp(1:i-1);
end

X = Q'*xp;


