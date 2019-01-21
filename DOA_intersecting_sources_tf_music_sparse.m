%computing MSE curve vs. SNR when 50% samples are randomly removed along with comparison with TF MUSIC
%In this code we assume that the user add in this path the TFSAP toolboox.
%Please cite the following reference:


close all;
clear all;
N_sensors=8;
n=0:127;
perc    = 0.1;

%addpath('D:\tfsa_5-5\windows\win64_bin');
%crossing componentsi8

s1=exp(2*pi*1i*(0.05*n+0.1*n.^2/(2*128)+0.2*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.1*n+0.1*n.^2/(2*128)+0.2*n.^3/(128*128*3)));
s3=exp(2*pi*1i*(0.4*n+0.1*n.^2/(2*128)-0.5*n.^3/(128*128*3)));%+exp(2*pi*1i*(0.45*n+0*0.1*n.^2/(2*128)-0*0.5*n.^3/(128*128*3)));
s4=exp(2*pi*1i*(0.45*n+0.1*n.^2/(2*128)-0.1*n.^3/(128*128*3)));
%crossing components LFM sources
s = [(s1.') (s2.') (s3.') (s4.')];
%s = [(s1.') (s2.') (s3.') ];

s1=exp(2*pi*1i*(0.05*n+0*0.2*n.^2/(2*128)));
s2=exp(2*pi*1i*(0.1*n+0.3*n.^2/(2*128)));
s3=exp(2*pi*1i*(0.4*n-0.3*n.^2/(2*128)));
s4=exp(2*pi*1i*(0.45*n-0.2*n.^2/(2*128)));

n_sources=4;
N_C=4;
s_orig=s;

% set mixing matrix A
theta = [15,30,50]*pi/180;   % sensor separation angles in radians
theta = [15,30,45,60]*pi/180;   % sensor separation angles in radians

A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A
X = A*s.';                             % mixed source
theta9=round(theta *180/pi);
% generate noise
SNR=10;
sigma = 10^(-SNR/20);
w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise

X=X+w;

for kk=1:N_sensors
    r= randsample(128,32*2) ;
    X(kk,r)=0;
end

% IF estimation
ss= multi_sensor_source_separation(X, N_C, 3,N_sensors);

%%DOA estimation
IP=1;
figure;
for iii=1:n_sources
    for jjj=1:N_sensors
        a(jjj,:)=ss(jjj,iii,:);
    end
    theta1=1:1:90;
    
    p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
    [x,y]=max(p);
    % y1(iii)=y;
    locs(iii)=y;
    P(iii,:)=p;
    % p=music((a), 1, N_sensors, 1, 1, theta1');
    % P(iii,:)=p;
end

plot(P')
aa=zeros(1,90*IP(1:IP:end));
aa(theta9)=-50;
aa(theta9)=1;
locs
hold on; stem(aa)
xlabel('DOA Estimation (degrees)');
ylabel('Spatial Spectrum');

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
[~,locs]=findpeaks(P_tf_music_ckd,'NPeaks',4,'MinPeakDistance',10,'Threshold',0.001);
if length(locs)<n_sources
    locs(locs+1:n_sources)=0;
elseif length(locs)>n_sources
    locs=locs(1:n_sources);
end
figure;plot(P_tf_music_ckd)
hold on; stem(aa)
xlabel('DOA Estimation (degrees)');
ylabel('Spatial Spectrum');
locs