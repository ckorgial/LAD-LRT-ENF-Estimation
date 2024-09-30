%% LAD ENF estimation based on search within sum of harmonic components %%%
% This function is based on the weighted maximum-likelihood estimator proposed
% by Hajj-Ahmad et al.
%
%   input:
%           signal:         time domain signal row vector
%           fs:             sampling frequency
%           window_dur:     window duration in second
%           step_size_dur:  step size duration in second
%           fc:             nominal multi-tone frequencies [100,150,200,...]
%           bound:          tolerable deviation of IF from fc in Hz
%           FFT_res_factor: number of FFT points = FFT_res_factor*fs
%           T:              sampling interval (1/fs)
%   output:
%           IF:            estimated ENF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IF,weights] = func_LAD_multi_tone_search_weighted(signal,fs,window_dur,step_size_dur,fc,bound,FFT_res_factor,T)
window_length   = window_dur*fs;
window_func     = rectwin(window_length)';
step_size       = step_size_dur*fs;
NFFT            = FFT_res_factor*fs;
window_pos      = 1:step_size:(length(signal)-window_length+1);
IF              = zeros(1,length(window_pos)); % output IF without interpolation
%% set bandwidth to estimate local SNR according to [2]
band_low_signal_1st   = round(49.98*FFT_res_factor);
band_high_signal_1st  = round(50.02*FFT_res_factor); % signal band: 49.98~50.02 Hz
band_low_noise_1st    = round(49.9*FFT_res_factor);
band_high_noise_1st   = round(50.1*FFT_res_factor); % noise band: 49.9~50.1 Hz excluding signal band                   
weights               = zeros(1,length(fc));
%% set harmonic search region
search_region_1st = round((50-bound(1)/2)*FFT_res_factor):round((50+bound(1)/2)*FFT_res_factor); % fundamental IF search region
length_per_band   = length(search_region_1st);
search_region     = kron(search_region_1st,(fc/50))'; % harmonic IF search region index
%% search loop
for i = 1:length(window_pos)
    temp           = fft(signal(window_pos(i):window_pos(i)+window_length-1).*window_func,NFFT);
    %abstemp        = signal(window_pos(i):window_pos(i)+window_length-1).*window_func;
    %N              = length(abstemp);
    N              = length(temp);
    HalfTempFFT    = temp(1:end/2);
    absHalfTempFFT = abs(HalfTempFFT).';
    % calculate weights for each harmonic component according to [2]
    for j = 1:length(fc)
       signal_component   = absHalfTempFFT(band_low_signal_1st*(fc(j)/50):band_high_signal_1st*(fc(j)/50));
       signal_plus_noise  = absHalfTempFFT(band_low_noise_1st*(fc(j)/50):band_high_noise_1st*(fc(j)/50));
       weights(j)         = norm(signal_component).^2/(norm(signal_plus_noise).^2-norm(signal_component).^2);
    end
    weights        = weights/norm(weights);
    fbin_candidate = absHalfTempFFT(search_region);
    fbin_candidate = diag(weights)*reshape(fbin_candidate,[length(fc),length_per_band]);
    weighted_fbin  = sum(fbin_candidate.^2,1); % harmonically weighted frequency bin energy 
    ValueMax       = max(weighted_fbin);
    PeakLoc        = search_region_1st((weighted_fbin==ValueMax(1)));
    % Non-Linear LAD
    convergence_threshold = 1e-4; % Set it properly
    iteration = 0;
    CONT=1;
    while CONT 
        iteration = iteration + 1;
        if iteration == 1
            fcc = PeakLoc;
        end
        % Optimization wrt theta
        Hm    = [cos(2*pi*T*fcc*(0:N-1))',sin(2*pi*T*fcc*(0:N-1))'];
        theta = ladreg(temp', Hm, false, [], 1); 
        % Optimization wrt fm
        zmin=1e7;
        for m = 1:99
            fm  = fcc + (m-50)*fs/(60*N);  
            Hm  = [cos(2*pi*T*fm*(0:N-1))',sin(2*pi*T*fm*(0:N-1))'];
            zm  = norm(temp' - Hm*theta,1); % Fix theta
            if zm < zmin
                fcc_new=fm;
                zmin=zm; 
            end   
        end
        relative_difference = abs(fcc - fcc_new)/fcc ;
        if relative_difference > convergence_threshold
            fcc=fcc_new;
        else
            IF(i)=fcc*fs/NFFT*2;
            CONT=0;
        end
    end
end
norm_fc            = 100;
norm_bound         = 100*bound(1)/fc(1);
IF(IF<norm_fc-norm_bound)=norm_fc-norm_bound;IF(IF>norm_fc+norm_bound)=norm_fc+norm_bound;
end