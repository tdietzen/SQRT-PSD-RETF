%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2019 Thoms Dietzen
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt).
%
% If you find it useful, please cite:
%
% [1] T. Dietzen, S. Doclo, M. Moonen, and T. van Waterschoot, “Square
% rootbased multi-source early PSD estimation and recursive RETF update in
% reverberant environments by means of the orthogonal Procrustes problem,”
% ESAT-STADIUS Tech. Rep. TR 19-69, KU Leuven, Belgium, submitted for
% publication, June 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of early PSD estimation as described in [1]. The code contained
% in main.m loads microphone signals and estimates early PSDs, based on
% both the conventional and the square-root MP. If the square-root MP is
% used, RETF updates are performed. The performance is measured in terms of
% SIR, SAR, and SDR.


%% PREAMBLE SETTINGS

clear;
cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
set(0,'DefaultFigureWindowStyle','docked');


%% CONFIGURATION

%%% ACOUSTIC SETTINGS

% speed of sound
c = 340;
% microphone positions
micPos = [...
    0, 0;...
    0.08, 0;...
    0.16, 0;...
    0.24, 0;...
    0.32, 0;...
    ];
% number of microphones
M = size(micPos,1);
% source angles
sourceAng = [0 60];
% microphone signals
x1_TD = audioread('x1.wav');
x2_TD = audioread('x2.wav');
s1_TD = audioread('s1.wav');
s2_TD = audioread('s2.wav');
y_TD = x1_TD + x2_TD;


%%% ALGORITHMIC SETTINGS

% sample rate
fs = 16000;
% STFT parameters
N_STFT = 512;
R_STFT = N_STFT/2;
win = sqrt(hann(N_STFT,'periodic'));
N_STFT_half = floor(N_STFT/2)+1;
% frequency vector
f = linspace(0,fs/2,N_STFT_half);

% forgetting factor zeta
zeta  = tau2forget(2*M*R_STFT/fs, R_STFT, fs);
% laplace coefficients of speech per frequency bin
tmp          = load('./source/param/lap_div.mat');
lap_div      = tmp.lap_div;
% penelty factor alpha for square-root MP
tmp          = load('./source/param/alpha_sqrtMP.mat');
alpha_sqrtMP = tmp.alpha_sqrtMP;

% yhird-octave band settings
numBands    = 18;
kappa       = -(numBands-1):1:0;
f_high      = 2.^(kappa/3)*f(end);
f_low       = 2.^((kappa-1)/3)*f(end);
k_low       = zeros(numBands,1);
k_high      = zeros(numBands,1);
for i_band = 1:numBands
    [~, k_low(i_band)]  = min(abs(f-f_low(i_band)));
    [~, k_high(i_band)] = min(abs(f-f_high(i_band)));
end

    
%%% INITIAL RETFs,  DIFFUSE COHERENCE MATRIX

% initial RETFs
H_init_FT = doa2steervec(micPos, sourceAng, N_STFT_half, fs, c);
% diffuse coherence matrix
Gamma = calc_diffcoherence(micPos,N_STFT,fs,c,1e-3);


%%% FIGURE SETTINGS

% third-octave band figure settings
TOxTicks = 0.5:3:18.5;
TOxTicklabels = 2.^(-numBands/3:1:0)*f(end)/1000;

% spectogram figure settings
xTickProp = [0, R_STFT/fs, fs/R_STFT];
yTickProp = [0, fs/(2000*R_STFT), R_STFT/2];
cRange    = [-55 5];


%% STFT PROCESSING

% transform
x1_STFT = calc_STFT(x1_TD, fs, win, N_STFT, R_STFT, 'onesided');
x2_STFT = calc_STFT(x2_TD, fs, win, N_STFT, R_STFT, 'onesided');
s1_STFT = calc_STFT(s1_TD, fs, win, N_STFT, R_STFT, 'onesided');
s2_STFT = calc_STFT(s2_TD, fs, win, N_STFT, R_STFT, 'onesided');
y_STFT  = x1_STFT + x2_STFT;
% plot
figure('Name','microphone signals');
subplot(3,1,1); plotSpec(x1_STFT(:,:,1),  'mag', [],        yTickProp, cRange, 0); title('x1'); ylabel('f/kHz');
subplot(3,1,2); plotSpec(x2_STFT(:,:,1),  'mag', [],        yTickProp, cRange, 0); title('x2'); ylabel('f/kHz');
subplot(3,1,3); plotSpec(y_STFT(:,:,1),   'mag', xTickProp, yTickProp, cRange, 0); title('y');  ylabel('f/kHz'); xlabel('time/s');
drawnow;


%% SUBSPACE PROCESSING

% correlation matrix of microphone signal
Psi_x_STFT = estim_corrmat(y_STFT, zeta);
% compute GEVD
[P_STFT, lambda_STFT] = desmooth_GEVD(Psi_x_STFT, Gamma,...
    'lambdaMin', 0,...
    'forgetPSD', zeta);


%% ESTIMATE PSDs, RETFs

%%% INIT

% estimates
phi_s_hat           = cell(3,1);
phi_xl_hat          = cell(3,1);
H_hat_prior_STFT    = cell(3,1);
H_hat_post_STFT     = cell(3,1);
H_update_pattern    = cell(3,1);

% measures
sqrt_phi_s_bar      = cell(2,1);
e_phi_s_int         = cell(2,1);
e_phi_s_art         = cell(2,1);
measures            = cell(2,1);

% reference PSD
phi_s_ref_STFT = cat(3, abs(s1_STFT).^2, abs(s2_STFT).^2);

for i_method = 1:2
    
    %%% SETTINGS
    
    switch i_method
        case 1
            method = 'conventional MP';
            alpha = zeros(N_STFT_half,1);
            beta = [];
            xi_thresh = [];
        case 2
            method = 'square-root MP';
            alpha = alpha_sqrtMP;
            beta = 20*lap_div.^2;
            xi_thresh = db2pow(-2);
    end
    
    
    %%% ESTIMATE PSDs, UPDATE RETFS
    
    [phi_s_hat{i_method},...
        phi_xl_hat{i_method},...
        H_hat_prior_STFT{i_method},...
        H_hat_post_STFT{i_method},...
        H_update_pattern{i_method}]...
        = estim_PSD_RETF(P_STFT, lambda_STFT, Gamma, H_init_FT,...
        'method', method,...
        'itmax', 20,...
        'alpha', alpha,...
        'beta', beta,...
        'xiThresh', xi_thresh...
        );
    

    %%% EVALUATE
    
    [sqrt_phi_s_bar{i_method},...
        e_phi_s_int{i_method},...
        e_phi_s_art{i_method},...
        measures{i_method},...
        ] = bss_eval_sources_wrapper(...
        sqrt(phi_s_hat{i_method}),...
        sqrt(phi_s_ref_STFT),...
        k_low,...
        k_high,...
        1);
    
    
    %%% PLOT
    
    figure('Name', [method ' PSD estimate - target']);
    subplot(3,1,1); plotSpec(sqrt_phi_s_bar{i_method}(:,:,1).^2,   'pow', [],        yTickProp, cRange, 0); title('phi s1 bar');   ylabel('f/kHz');
    subplot(3,1,2); plotSpec(sqrt_phi_s_bar{i_method}(:,:,2).^2,   'pow', [],        yTickProp, cRange, 0); title('phi s2 bar');   ylabel('f/kHz');
    subplot(3,1,3); plotSpec(sum(sqrt_phi_s_bar{i_method}.^2,3),   'pow', xTickProp, yTickProp, cRange, 0); title('phi tot bar'); ylabel('f/kHz');

    figure('Name', [method ' PSD estimate - interfecence']);
    subplot(3,1,1); plotSpec(e_phi_s_int{i_method}(:,:,1).^2,      'pow', [],        yTickProp, cRange, 0); title('e phi s1 int'); ylabel('f/kHz');
    subplot(3,1,2); plotSpec(e_phi_s_int{i_method}(:,:,2).^2,      'pow', [],        yTickProp, cRange, 0); title('e phi s2 int'); ylabel('f/kHz');
    subplot(3,1,3); plotSpec(sum(e_phi_s_int{i_method}.^2,3),      'pow', xTickProp, yTickProp, cRange, 0); title('e phi tot int'); ylabel('f/kHz');

    figure('Name', [method ' PSD estimate - artifacts']);
    subplot(3,1,1); plotSpec(e_phi_s_art{i_method}(:,:,1).^2,      'pow', [],        yTickProp, cRange, 0); title('e phi s1 art'); ylabel('f/kHz');
    subplot(3,1,2); plotSpec(e_phi_s_art{i_method}(:,:,2).^2,      'pow', [],        yTickProp, cRange, 0); title('e phi s2 art'); ylabel('f/kHz');
    subplot(3,1,3); plotSpec(sum(e_phi_s_art{i_method}.^2,3),      'pow', xTickProp, yTickProp, cRange, 0); title('e phi tot art'); ylabel('f/kHz');

    if strcmp(method, 'square-root MP')
        figure('Name', [method ' RETF update pattern']);
        subplot(3,1,1); plotSpec(H_update_pattern{i_method}(:,:,1),     'lin', [],        yTickProp, [0, 1], 0);   title('update h1'); ylabel('f/kHz');
        subplot(3,1,2); plotSpec(H_update_pattern{i_method}(:,:,2),     'lin', xTickProp, yTickProp, [0, 1], 0);   title('update h2'); ylabel('f/kHz');
        drawnow;
    end
    
    figure('Name', [method ' performance measures']);
    subplot(3,1,1); bar(measures{i_method}.SIR_band); xlim([0 19]); ylim([0 22]); ylabel('SIR/dB'); set(gca,'xTick',[]);
    subplot(3,1,2); bar(measures{i_method}.SAR_band); xlim([0 19]); ylim([0 22]); ylabel('SAR/dB'); set(gca,'xTick',[]);
    subplot(3,1,3); bar(measures{i_method}.SDR_band); xlim([0 19]); ylim([0 22]); ylabel('SDR/dB'); set(gca,'xTick',TOxTicks); set(gca,'xTickLabel',TOxTicklabels); xlabel('f/kHz');
    drawnow;
    
end
