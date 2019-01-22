%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] B. Boashash, A. Aissa-El-Bey, Multisensor Time-Frequency Signal Processing:
%     A tutorial review with illustrations in selected application areas, Digital
%     Signal Processing, In Press.
% [2] B. Boashash, A. Aissa-El-Bey, M. F. Al-Sa'd, Multisensor time-frequency
%     signal processing software Matlab package: An analysis tool for multichannel
%     non-stationary data , SoftwareX, In Press.
%
% Last Modification: 21-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Short Time-Fourier Transform (STFT)
%
% Syntax : TFD_stft = stft(s, M, h);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
% <INPUTs>
% s    : Input signal. This can be real or analytic.
% M    : Number of Frequency bins (FFT Length).
% h    : The Lag window name. See help of the Matlab function "window" to
%        list all available windows.
%
% <OUTPUTs>
% TFD_stft : The STFT of the input signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% s1 = chirp(0:255,0.4,255,0.1);
% s2 = chirp(0:255,0.1,255,0.4);
% s = s1 + s2;
% TFD = stft(s, 512, hamming(63)');
% figure; imagesc(0:1/1023:1/2,0:255,abs(TFD')); axis xy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TFD_stft = stft(s, M, h)
%% Check for Errors
[L1, M, N] = input_logical_errors(s, length(h), M);
L = (L1-1)/2;

%% Computing the Signal Kernel and Smoothing in the Time-Lag Domain
K_TL = zeros (M, N);
for n = 1:N,
    tau = -min([round(M/2)-1,L,n-1]):min([round(M/2)-1,L,N-n]);
    ind = rem(M+tau,M)+1;
    K_TL(ind,n)=s(n+tau).*conj(h(L+1+tau))/norm(h(L+1+tau));
end
TFD_stft = fft(K_TL);
TFD_stft = TFD_stft(1:end/2,:);

    %% Check for Logical Errors
    function [L, M, N] = input_logical_errors(s, L, M)
        N = length(s);
        if L >= N
            L = N - 1; disp('Window length is truncated to N-1')
        end
        if M < N
            M = N; disp('FFT length is extended to N')
        end
        if rem(L,2) == 0
            L = L - 1; disp('Window length must be odd, thus it is truncated to L-1')
        end
        if(L <= 0)
            error('Window length must be positive');
        end
        if(M <= 0)
            error('FFT length must be positive');
        end
        M = 2*M;
    end
end