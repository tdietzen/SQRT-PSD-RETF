function [ target, interf, artif, measure ] = bss_eval_sources_wrapper(sqrt_phi_s_hat_STFT, sqrt_phi_s_STFT, k_low, k_high, flen)
% [ target, interf, artif, measure ] = bss_eval_sources_wrapper(sqrt_phi_s_hat_STFT, sqrt_phi_s_STFT, k_low, k_high, flen)
% wrapper for bss_eval_sources_edit.
%
% IN:
% sqrt_phi_s_hat_STFT       square root of early PSD estimates - freqbins x frames x sources
% sqrt_phi_s_STFT           square root of early PSD reference - freqbins x frames x sources
% k_low                     freqbin indices for lower end of third octave bands
% k_high                    freqbin indices for upper end of third octave bands
% flen                      filter length in bss_eval_sources_edit
%
% OUT:
% target                    PSD estimate target component - freqbins x frames x sources
% interf                    PSD estimate interference component - freqbins x frames x sources
% artif                     PSD estimate artifact component - freqbins x frames x sources
% measure                   struct containing SIR, SAR, and SDR scores

[N_FT_half, L, N] = size(sqrt_phi_s_STFT);

target = zeros(N_FT_half,L,N);
interf = zeros(N_FT_half,L,N);
artif  = zeros(N_FT_half,L,N);

for k = 1:N_FT_half

    sqrt_phi_s      = squeeze(sqrt_phi_s_STFT(k,:,:)).';
    sqrt_phi_s_hat  = squeeze(sqrt_phi_s_hat_STFT(k,:,:)).';
    
    [~, ~, ~, target_tmp, interf_tmp, artif_tmp] = bss_eval_sources_edit(sqrt_phi_s_hat, sqrt_phi_s, flen);
    target(k,:,:) = target_tmp.';
    interf(k,:,:) = interf_tmp.';
    artif(k,:,:)  = artif_tmp.';
    
end

% N_FT_half x N_src
pow_target          = squeeze(mean(abs(target).^2,2));
pow_interf          = squeeze(mean(abs(interf).^2,2));
pow_artif           = squeeze(mean(abs(artif).^2,2));

% N_bands x 1
numBands = length(k_low);
pow_target_band = zeros(numBands,1);
pow_interf_band = zeros(numBands,1);
pow_artif_band  = zeros(numBands,1);
for i_band = 1:numBands
    if i_band == numBands
       i_bandrange = k_low(i_band):k_high(i_band);
    else
       i_bandrange = k_low(i_band):k_high(i_band)-1;
    end
    pow_target_band(i_band) = sum(sum(pow_target(i_bandrange,:)));
    pow_interf_band(i_band) = sum(sum(pow_interf(i_bandrange,:)));
    pow_artif_band(i_band)  = sum(sum(pow_artif(i_bandrange,:)));
end
% N_FT_half x 1
SDR_band  = pow2db(pow_target_band./                      (pow_interf_band+pow_artif_band    +1e-12));
SIR_band  = pow2db(pow_target_band./                      (pow_interf_band                   +1e-12));
SAR_band  = pow2db((pow_target_band+pow_interf_band)./    (pow_artif_band                    +1e-12));

%%%%%%% SAVE %%%%%%%%%%%%

% N_FT_half x N_src
measure.pow_target = pow_target;         
measure.pow_interf = pow_interf;       
measure.pow_artif  = pow_artif;          

% N_bands x 1
measure.pow_target_band = pow_target_band;         
measure.pow_interf_band = pow_interf_band;       
measure.pow_artif_band  = pow_artif_band;  

% N_bands x 1
measure.SDR_band = SDR_band; 
measure.SIR_band = SIR_band; 
measure.SAR_band = SAR_band;


end