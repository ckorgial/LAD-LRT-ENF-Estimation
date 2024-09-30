%% LAD single-tone ENF estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   input:
%           signal:         time domain signal row vector
%           fs:             sampling frequency
%           window_dur:     window duration in second
%           step_size_dur:  step size duration in second
%           fc:             nominal single tone frequency (50 Hz, 100 Hz, etc.)
%           bound:          tolerable deviation of IF from fc in Hz
%           FFT_res_factor: number of FFT points = FFT_res_factor*fs
%           T:              sampling interval (1/fs)
%   output:
%           IF:             estimated ENF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IF]   = func_LAD_single_tone(signal,fs,window_dur,step_size_dur,fc,bound,FFT_res_factor,T)
window_length   = window_dur*fs; % 16*fs
window_func     = rectwin(window_length)';
step_size       = step_size_dur*fs; % 1*fs
NFFT            = FFT_res_factor*fs; 
window_pos      = 1:step_size:(length(signal)-window_length+1); % length(input) = 386165
IF              = zeros(1,length(window_pos)); % output IF without interpolation
absHalfTempFFT0 = zeros(1,NFFT/2);
N_in1           = round(NFFT/fs);
N_in2           = round(0.25*NFFT/fs);
for i = 1:length(window_pos)
    temp            = fft(signal(window_pos(i):window_pos(i)+window_length-1).*window_func,NFFT);
    N               = length(temp);
    HalfTempFFT     = temp(1:end/2);
    absHalfTempFFT  = abs(HalfTempFFT);
    absHalfTempFFT0(fc*N_in1-N_in2:fc*N_in1+N_in2) = absHalfTempFFT(fc*N_in1-N_in2:fc*N_in1+N_in2);
    ValueMax        = max(absHalfTempFFT0);
    PeakLoc         = find(absHalfTempFFT0 == ValueMax(1))*(fs/NFFT);
    PeakLoc         = PeakLoc(1);   
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
            IF(i)=fcc;
            CONT=0;
        end
    end
end
IF(IF<fc-bound)=fc-bound;IF(IF>fc+bound)=fc+bound;
end