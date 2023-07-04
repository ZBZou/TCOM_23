function hat_r = division(rx,M)

K = log2(M);
Nu = size(rx,1);
min = -1-2*(K-1);
max = 1+2*(K-1);
hat_r = rx;
for i = 1:Nu
   hat_r(i,:) = round(rx(i,:));
   hat_r(i,find(hat_r(i,:) < min)) = min;
   hat_r(i,find(hat_r(i,:) > max)) = max;
end

end