%computing MSE curve in over-determined scenario along with comparison with TF MUSIC
%In this code we assume that the user add in this path the TFSAP toolboox.
%Please cite the following reference:

close all;
clear all;

% NUmber of simulation runs
LL=200;

% Signal model%
addpath 'D:\tfsa_5-5\windows\win64_bin\'
index=0;
n=0:127;
% Number of SOurces
n_sources=4;
% NUmber of components
N_C=n_sources;
N_sensors=8;
s1=exp(2*pi*1i*(0.05*n+0.1*n.^2/(2*128)+0.2*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.1*n+0.1*n.^2/(2*128)+0.2*n.^3/(128*128*3)));
s3=exp(2*pi*1i*(0.4*n+0.1*n.^2/(2*128)-0.5*n.^3/(128*128*3)));%+exp(2*pi*1i*(0.45*n+0*0.1*n.^2/(2*128)-0*0.5*n.^3/(128*128*3)));
s4=exp(2*pi*1i*(0.45*n+0.1*n.^2/(2*128)-0.1*n.^3/(128*128*3)));

perc=0.4;
% s1=exp(2*pi*1i*(0.05*n+0*0.2*n.^2/(2*128)));
% s2=exp(2*pi*1i*(0.1*n+0.3*n.^2/(2*128)));
% s3=exp(2*pi*1i*(0.4*n-0.3*n.^2/(2*128)));
% s4=exp(2*pi*1i*(0.45*n-0.2*n.^2/(2*128)));

for SNR=-5:1:10
    
    for ii=1:LL
        M=N_sensors;% Number of sensors
        
        
        %Close components
        
        %s4 =exp(2*pi*1i*(0.35*n-1.25*1.0863e-06*n.^3));
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
        %% The Proposed algorithm
        %Separation of sources using TF filtering
        
        
        
        
        ss= multi_sensor_source_separation(X, N_C, 2,N_sensors);
        
        
        
        
        
        %K=K1;
        %
        
        
        
        
        
        
        
        
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
        
        mmssee_adtfd(ii)=mean((sort(y1)-sort(theta)).^2);
        
        % The proposed algorithm ends
        
        D   = mtfd(X, 'ckd',1, 0.3, 0.3, length(X));
        %%% Averaged Auto-TFD
        D_avg = zeros(length(X), length(X));
        for mm = 1:N_sensors, D_avg = D{mm,mm} + D_avg; end
        D_avg = D_avg./N_sensors;
        %%% Selection of high-energy (t,f) points
        thr = perc*max(max(D_avg));
        Tr = abs(D_avg) >= thr;
        [F_trace, ~] = find(Tr);
        n_p = length(F_trace);
        D_s = zeros(N_sensors, N_sensors);
        for m1 = 1:N_sensors
            for m2 = 1:N_sensors
                D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
            end
        end
        %%% DOA Estimation
        %P_tf_music_ckd = tf_music(D_s, n_sources, N_sensors, 2,1, theta1);
        P_tf_music_ckd =TMMUSIC(D_s, 2, N_sensors, n_sources, 1, theta1);
        [~,y1]=findpeaks(P_tf_music_ckd,'NPeaks',4,'MinPeakDistance',10,'Threshold',0.001);
        %figure;plot(P_tf_music_ckd)
        if length(y1)<n_sources
            y1(length(y1)+1:n_sources)=0;
        elseif length(y1)>n_sources
            y1=y1(1:n_sources);
        end
        y1=3*y1/180;
        
        mmssee_tf_music(ii)=mean((sort(y1)-sort(theta)).^2);
        
        
    end
    index=index+1;
    mean(mmssee_adtfd)
    mean(mmssee_tf_music)
    
    snr_mse_adtfd(index)=mean(mmssee_adtfd);
    snr_mse_tf(index)=mean(mmssee_tf_music);
end
snr_mse_adtfd(:)
snr_mse_tf(:)
index=0;



SNR=-5:1:10;
plot(SNR(1:5:end),10*(log10(snr_mse_adtfd(1:5:end))),'--md','linewidth',2);
hold on;
plot(SNR(1:5:end),10*(log10(snr_mse_tf(1:5:end))),'r','linewidth',2);
xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
legend('The Proposed Algorithm','Time-frequency Music');
% axis([min(snr) max(snr)  -50  0])

%    snr_mse(1,5)=snr_mse(1,4);
%
%     figure;
%     plot(SNR(1:end),1*((snr_mse(1,1:end))),'-ro','linewidth',2);
% %
%     hold on;
%     plot(SNR(1:end),1*((snr_mse(2,1:end))),':gs','linewidth',2);
%
%     hold on;
%     plot(SNR(1:end),1*((snr_mse(3,1:end))),'-.b+','linewidth',2);
%    hold on;
%    plot(SNR(3:end),1*((snr_mse(4,3:end))),'--md','linewidth',2);
%
%%    xlabel('Signal to Noise Ratio');
%    ylabel('Mean Square Error (dB)');



