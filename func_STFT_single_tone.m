%% Single-tone ENF estimation. Values of IFs bounded by [fc-bound, fc+bound]
%% Use this function only if input signal is properly bandpass filtered.
%   input:
%           signal:         time domain signal row vector
%           fs:             sampling frequency
%           window_dur:     window duration in second
%           step_size_dur:  step size duration in second
%           fc:             nominal single tone frequency (50 Hz, 100 Hz, etc.)
%           bound:          tolerable deviation of IF from fc in Hz
%           FFT_res_factor: number of FFT points = FFT_res_factor*fs
%   output:
%           IF:            estimated ENF wihtout interpolation
%           STFT_TFD:      estimated time-frequency distribution (TFD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IF,STFT_TFD] = func_STFT_single_tone(signal,fs,window_dur,step_size_dur,fc,bound,FFT_res_factor)
window_length   = window_dur*fs;
window_func     = rectwin(window_length)';
step_size       = step_size_dur*fs;
NFFT            = FFT_res_factor*fs;
window_pos      = 1:step_size:(length(signal)-window_length+1);
STFT_TFD        = zeros (NFFT,length(window_pos));
IF              = zeros(1,length(window_pos)); % output IF without interpolation
absHalfTempFFT0 = zeros(1,NFFT/2);
N_in1           = round(NFFT/fs);
N_in2           = round(0.25*NFFT/fs);
for i = 1:length(window_pos)
    temp            = fft(signal(window_pos(i):window_pos(i)+window_length-1).*window_func,NFFT);
    STFT_TFD(:,i)   = temp;
    HalfTempFFT     = temp(1:end/2);
    absHalfTempFFT  = abs(HalfTempFFT);
    absHalfTempFFT0(fc*N_in1-N_in2:fc*N_in1+N_in2) = absHalfTempFFT(fc*N_in1-N_in2:fc*N_in1+N_in2);
    ValueMax        = max(absHalfTempFFT0);
    PeakLoc         = find(absHalfTempFFT0==ValueMax(1));
    PeakLoc         = PeakLoc(1);
    IF(i)           = (PeakLoc-1)*fs/NFFT;
end
IF(IF<fc-bound)=fc-bound;IF(IF>fc+bound)=fc+bound;
end