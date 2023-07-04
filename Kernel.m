%%%%%%kH%%%%
function kH = Kernel(LH, T1, F1, T2, F2)
kH = zeros(T1, F1, T2, F2);
for i = 1:T1
   for j = 1:F1
      for n = 1:T2
         for m = 1:F2
            if i-n < 1 || j-m < 1
               continue
            end
            kH(i,j,n,m) = LH(i-n,j-m)*exp(-2i*pi*(j - m)*n);
         end
      end
   end
end 
