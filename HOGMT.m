function X_hogmt = HOGMT(s, h, ep)
Nep = length(ep);
Nu = size(h,1);
T = size(h,2);
X_hogmt = zeros(Nu, T, Nep);
H = reshape(h, Nu*T, Nu*T);
S = reshape(s,Nu*T,1);
display("Decomposing Channel");
[U,temp_sig,V] = svd(H);
display("Decomposed");
for i = 1:Nep
    Ne = round(Nu*T*ep(i));
    X = V(:,1:Ne)*pinv(temp_sig(1:Ne,1:Ne))*U(:,1:1:Ne)'*S;
    % X = V*pinv(temp_sig)*U'*S;
    X_hogmt(:,:,i) = reshape(X,Nu,T);
end

% Sig = diag(temp_sig);
% for i = 1:length(Sig)                
%    if Sig(i) > ep
%       N = i; % N most contributed eigenfunctions
%    end
% end 
% 
% % construct X
% S = reshape(s,Nu*T,1);
% 
% X = zeros(length(S),1);
% 
% for i = 1:N
%    xn = dot(S, U(:,i))/Sig(i);
%    X = X + xn*(V(:,i));
% end
% 
% X = reshape(X,Nu,T);

% Nu = size(h,1);
% T = size(h,2);
% h = reshape(h, Nu*T, Nu*T);
% s = reshape(s,[],1);
% h_r = real(h);
% h_i = imag(h);
% H = [h_r, -h_i;h_i, h_r];
% display("Decomposing Channel");
% [U,temp_sig,V] = svd(H);
% display("Decomposed");
% Sig = diag(temp_sig);
% for i = 1:length(Sig)                
%    if Sig(i) >= ep
%       N = i; % N most contributed eigenfunctions
%    end
% end 
% 
% % construct X
% sr = real(s);
% si = imag(s);
% 
% S = [sr;si];
% 
% X = zeros(length(S),1);
% 
% for i = 1:N
%    xn = dot(S, U(:,i))/Sig(i);
%    X = X + xn*V(:,i);
% end
% 
% Xr = X(1:length(s));
% Xi = X(length(s)+1:end)*1i;
% X = Xr+Xi;
% X = reshape(X,Nu,T);
