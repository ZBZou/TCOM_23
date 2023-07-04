function [ecdf,pcdf,scdf] = eigencdfs(s, h)
Nu = size(h,1);
T = size(h,2);
H = reshape(h, Nu*T, Nu*T);
S = reshape(s,Nu*T,1);
X = zeros(Nu*T,1);
display("Decomposing Channel");
[U,temp_sig,V] = svd(H);
display("Decomposed");
%%
sig = diag(temp_sig);
Nsig = length(sig);
ecdf = zeros(Nsig,1);
pcdf = zeros(Nsig,1);
scdf = zeros(Nsig,1);
s_coeff = conj(S')*U;
p = abs(s_coeff./conj(sig')).^2;
for i = 1:Nsig
   ecdf(i,1) = sum(sig(1:i));
   scdf(i,1) = sum(abs(s_coeff(1:i)));
   pcdf(i,1) = sum(p(1:i));
end
ecdf = ecdf/sum(sig);
pcdf = pcdf/sum(p);
scdf = scdf/sum(abs(s_coeff));
%%
plot([1:Nsig]/8000,ecdf,[1:Nsig]/8000,scdf,[1:Nsig]/8000,pcdf);
xlabel('Eigenfunctions');
ylabel('CDF');
legend('Eigenvalue','Reconstruction','Power');
