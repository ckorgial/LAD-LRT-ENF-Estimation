%% LAD ENF harmonic estimation algorithms with harmonic enhancement %%%%%%%
%  This MATLAB programm presents the LAD single/multi-tone harmonic model 
%  based ENF estimation algorithms.
% The code structure and implementation details are inspired by the
% methodology presented by Hua et al.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
tic;
%% parameter setting and initialization
FS                       = 800; % constant sampling frequency
T                        = 1/FS;
HARMONIC_INDEX           = [2,3,4,5,6,7]; % constant value for ENF harmonic processing
fc                       = 50*HARMONIC_INDEX; % nominal frequencies at each harmonic
bound                    = 0.1*HARMONIC_INDEX; % tolerable IF deviations at each harmonic

filter_length            = 256;
[BPF_coeffs, coeffs_2nd] = func_BPF(filter_length);
index                    = dir('C:\Users\30694\Documents\PHD\Codes\ENF_Enhancement_Estimation\Datasets\Recordings\H1\*.wav');
index_ref                = dir('C:\Users\30694\Documents\PHD\Codes\ENF_Enhancement_Estimation\Datasets\Recordings\H1_ref\*.wav');

load('selected_harmonic_index.mat')
load('selected_harmonic_index0.mat')

bound_MLE                = selected_harmonic_index; % cell array - output of the GHSA
bound_MLE                = cellfun(@(x) x*0.1, bound_MLE, 'UniformOutput', false); % tolerable IF deviations at each harmonic

bound_WMLE               = selected_harmonic_index0; % cell array - output of the GHSA
bound_WMLE               = cellfun(@(x) x*0.1, bound_WMLE, 'UniformOutput', false); % tolerable IF deviations at each harmonic

f_ref                    = cell(1,length(index));
 
f_LAD_single             = cell(1,length(index));
f_LAD_E_single           = cell(1,length(index));
f_LAD_MLE                = cell(1,length(index));
f_LAD_WMLE               = cell(1,length(index));
f_LAD_E_MLE              = cell(1,length(index));
f_LAD_E_WMLE             = cell(1,length(index));
f_LAD_S_MLE              = cell(1,length(index));
f_LAD_S_WMLE             = cell(1,length(index));
f_LAD_P_MLE              = cell(1,length(index));
f_LAD_P_WMLE             = cell(1,length(index));

laplacian_LAD_P_MLE    = cell(1,length(index));

MSE_LAD_single           = zeros(1,length(index));
MSE_LAD_E_single         = zeros(1,length(index));
MSE_LAD_MLE              = zeros(1,length(index));
MSE_LAD_WMLE             = zeros(1,length(index));
MSE_LAD_E_MLE            = zeros(1,length(index));
MSE_LAD_E_WMLE           = zeros(1,length(index));
MSE_LAD_S_MLE            = zeros(1,length(index));
MSE_LAD_S_WMLE           = zeros(1,length(index));
MSE_LAD_P_MLE            = zeros(1,length(index));
MSE_LAD_P_WMLE           = zeros(1,length(index));

tic;
for i=1:length(index)
    disp(['i=',num2str(i)]);
    %% import audio recording and reference  
    [audio, fs_audio] = audioread(strcat('C:\Users\30694\Documents\PHD\Codes\ENF_Enhancement_Estimation\Datasets\Recordings\H1\',index(i).name));
    [ref, fs_ref]     = audioread(strcat('C:\Users\30694\Documents\PHD\Codes\ENF_Enhancement_Estimation\Datasets\Recordings\H1_ref\',index_ref(i).name));
    ref               = ref';
    audio             = audio(:,1);
    audio             = audio';
    raw_wave          = resample(audio, FS, fs_audio);
    N                 = length(raw_wave);
    %% bandpass filtering
    input             = filtfilt(BPF_coeffs,1,raw_wave); % multi-tone input signal without enhancement
    %% harmonic enhancement
    N_ite             = 2; % enhancement iterations
    h_rfa             = 3000; % window length of RFA
    initial_guess     = fc(1)*ones(1,N); % initial IF of RFA is fixed to 100 Hz
    TS                = 1; % constant time-step for RFA, 1 second
    window_dur_rfa    = 8; % enhancement window size 8 seconds
    FFT_res_rfa       = 200; % FFT resolution for RFA 1/FFT_res_rfa Hz      
    refined_guess     = initial_guess;
    % step 1: single-tone enhancement
    for k = 1:N_ite
        [input_denoised_single,~,refined_guess] = func_RFA(input,h_rfa,FS,TS,...
            refined_guess,fc(1),bound(1),window_dur_rfa,FFT_res_rfa);
    end
    % step 2: use refined_guess to construct initial guesses for harmonics
    initial_guesses   = kron(HARMONIC_INDEX(2:end)'/2,refined_guess);
    refined_guesses   = initial_guesses;
    % step 3: harmonic enhancement
    for k=1:N_ite
        [input_denoised_multi,~,refined_guesses] = func_RFA_multi(input,h_rfa,FS,TS,...
            refined_guesses,fc(2:end),bound(2:end),window_dur_rfa,FFT_res_rfa);
    end
    input_E_single    = input_denoised_single;
    input_E_multi     = sum([input_denoised_single;input_denoised_multi],1);
    %% LAD ENF Estimators 
    % Set up parameters forframe-based processing
    window_dur        = 16; % duration of overlapping frame in second
    step_size_dur     = 1; % frame step-size usually 1 second
    FFT_res_factor    = 2000; % FFT resolution = 1/FFT_res_factor Hz
    f_ref{i}          = func_STFT_single_tone(ref,fs_ref,window_dur,step_size_dur,50,0.1,FFT_res_factor)*2; % reference
    % 1. single-tone estimation (2nd harmonic)
    f_LAD_single{i}   = func_LAD_single_tone(input,FS,window_dur,step_size_dur,fc(1),bound(1),FFT_res_factor,T);
    % 2. single-tone enhancement (2nd harmonic)
    f_LAD_E_single{i} = func_LAD_single_tone(input_E_single,FS,window_dur,step_size_dur,fc(1),bound(1),FFT_res_factor,T);
    % 3. Search within sum of harmonic components, mapped to 2nd harmonic
    f_LAD_MLE{i}      = func_LAD_multi_tone_search(input,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor,T);
    % 4. Search within weighted sum of harmonic components, mapped to 2nd harmonic
    f_LAD_WMLE{i}     = func_LAD_multi_tone_search_weighted(input,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor,T);
    % 5. Only enhancement of LAD_MLE
    f_LAD_E_MLE{i}    = func_LAD_multi_tone_search(input_E_multi,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor,T);
    % 6. Only enhancement of LAD_WMLE
    f_LAD_E_WMLE{i}   = func_LAD_multi_tone_search_weighted(input_E_multi,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor,T);
    % 7. LAD-S-MLE
    f_LAD_S_MLE{i}    = func_LAD_multi_tone_search(input,FS,window_dur,step_size_dur,fc,bound_MLE{i},2*FFT_res_factor,T);
    % 8. LAD-S-WMLE
    f_LAD_S_WMLE{i}   = func_LAD_multi_tone_search_weighted(input,FS,window_dur,step_size_dur,fc,bound_WMLE{i},2*FFT_res_factor,T);
    % 9. LAD-P-MLE
    [f_LAD_P_MLE{i}]  = func_LAD_multi_tone_search(input_E_multi,FS,window_dur,step_size_dur,fc,bound_MLE{i},2*FFT_res_factor,T);
    % 10. LAD-P-WMLE
    [f_LAD_P_WMLE{i}] = func_LAD_multi_tone_search_weighted(input_E_multi,FS,window_dur,step_size_dur,fc,bound_WMLE{i},2*FFT_res_factor,T);
    
    MSE_LAD_single(i)       = 1/length(f_ref{i})*norm(f_LAD_single{i}-f_ref{i}).^2;
    MSE_LAD_E_single(i)     = 1/length(f_ref{i})*norm(f_LAD_E_single{i}-f_ref{i}).^2;   
    MSE_LAD_MLE(i)          = 1/length(f_ref{i})*norm(f_LAD_MLE{i}-f_ref{i}).^2;
    MSE_LAD_WMLE(i)         = 1/length(f_ref{i})*norm(f_LAD_WMLE{i}-f_ref{i}).^2;
    MSE_LAD_E_MLE(i)        = 1/length(f_ref{i})*norm(f_LAD_E_MLE{i}-f_ref{i}).^2;
    MSE_LAD_E_WMLE(i)       = 1/length(f_ref{i})*norm(f_LAD_E_WMLE{i}-f_ref{i}).^2;
    MSE_LAD_S_MLE(i)        = 1/length(f_ref{i})*norm(f_LAD_S_MLE{i}-f_ref{i}).^2;
    MSE_LAD_S_WMLE(i)       = 1/length(f_ref{i})*norm(f_LAD_S_WMLE{i}-f_ref{i}).^2;
    MSE_LAD_P_MLE(i)        = 1/length(f_ref{i})*norm(f_LAD_P_MLE{i}-f_ref{i}).^2;
    MSE_LAD_P_WMLE(i)       = 1/length(f_ref{i})*norm(f_LAD_P_WMLE{i}-f_ref{i}).^2;
end
toc;

mean(MSE_LAD_single)
mean(MSE_LAD_E_single)
mean(MSE_LAD_MLE)
mean(MSE_LAD_WMLE)
mean(MSE_LAD_E_MLE)
mean(MSE_LAD_E_WMLE)

std(MSE_LAD_single)
std(MSE_LAD_E_single)
std(MSE_LAD_MLE)
std(MSE_LAD_WMLE)
std(MSE_LAD_E_MLE)
std(MSE_LAD_E_WMLE)