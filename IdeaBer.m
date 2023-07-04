function Ber_ideal = IdeaBer(data, s_ofdm, M, Nu, T, sigma, N_carrier, N_cp)
sigPow = sum(sum(abs(s_ofdm)))/(size(s_ofdm,1)*size(s_ofdm,2));
Rx_ideal = s_ofdm + sigPow*sigma*(randn(Nu,size(s_ofdm,2))+randn(Nu,size(s_ofdm,2))*1i)/sqrt(2);
r_temp = permute(Rx_ideal, [2,1]);
s_demod = ofdmdemod(r_temp, N_carrier, N_cp);
s_hat = reshape(permute(s_demod, [3,1,2]),Nu,[]);
Ber_ideal = myber(data, s_hat, M);