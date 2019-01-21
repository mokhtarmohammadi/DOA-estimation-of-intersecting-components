%% computing MSE curve vs. SNR in under-determined scenario
%In this code we assume that the user add in this path the TFSAP toolboox.
%Please cite the following reference:

close all;
clear all;
% NUmber of simulation runs
LL=100;

% Signal model%
%addpath 'D:\tfsa_5-5\windows\win64_bin\'
index=0;
n=0:127;
% Number of SOurces
n_sources=4;
% NUmber of components
N_C=n_sources;
N_sensors=3;
s1=exp(2*pi*1i*(0.05*n+0.1*n.^2/(2*128)+0.2*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.1*n+0.1*n.^2/(2*128)+0.2*n.^3/(128*128*3)));
s3=exp(2*pi*1i*(0.4*n+0.1*n.^2/(2*128)-0.5*n.^3/(128*128*3)));%+exp(2*pi*1i*(0.45*n+0*0.1*n.^2/(2*128)-0*0.5*n.^3/(128*128*3)));
s4=exp(2*pi*1i*(0.45*n+0.1*n.^2/(2*128)-0.1*n.^3/(128*128*3)));

for SNR=-5:1:10
    
    for ii=1:LL
        M=N_sensors;% Number of sensors
        
        % Signal generation
        s = [(s1.') (s2.') (s3.') (s4.')];
        s_orig=s;
        % set mixing matrix A
        theta = [10,30,50,70]*pi/180;   % sensor separation angles in radians
        A = exp(1j*pi*[0:M-1].'*sin(theta));  % mixing matrix A
        X = A*s.';                             % mixed source
        theta9=round(theta *180/pi);
        % generate noise
        sigma = 10^(-SNR/20);
        w = sigma*(randn(M,length(n)) + 1j*(randn(M,length(n))))/sqrt(2); % noise
        
        X=X+w;
        ss= multi_sensor_source_separation(X, N_C, 2,N_sensors);
        %%%%%%%%  BSS code  ends
        %%DOA estimation
        IP=1;
        %figure;
        clear y1;
        for iii=1:n_sources
            for jjj=1:N_sensors
                a(jjj,:)=ss(jjj,iii,:);
            end
            theta1=1:1:90;
            
            p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
            [x,y]=max(p);
             y1(iii)=y;
            P(iii,:)=p;
            % p=music((a), 1, N_sensors, 1, 1, theta1');
            % P(iii,:)=p;
        end
        
        if length(y1)>length(theta)
            y1=y1(1:4);
        elseif length(y1)<length(theta)
            y1(length(y1):length(theta))=0;
        end
        theta=3*theta/pi;
        y1=3*y1/180;
        
        mmssee(ii)=mean((sort(y1)-sort(theta)).^2);
    end
    index=index+1;
    mean(mmssee)
    snr_mse(index)=mean(mmssee);
end
snr_mse(:)
index=0;
SNR=-5:1:10;

plot(SNR(1:2:end),10*(log10(snr_mse(1:2:end))),'--md','linewidth',2);
xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
