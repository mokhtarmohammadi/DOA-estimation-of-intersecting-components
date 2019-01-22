%In this code we assume that the user add in this path the TFSAP toolboox.
%Please cite the following reference:

close all;
clear all;
N_sensors=3;
n=0:127;

%addpath('D:\D\win64_bin\win64_bin');
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
SNR=5;
sigma = 10^(-SNR/20);
w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise

X=X+w;

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
    P(iii,:)=p; 
  % p=music((a), 1, N_sensors, 1, 1, theta1');
   % P(iii,:)=p;
end

   plot(P')
aa=zeros(1,90*IP(1:IP:end));
aa(theta9)=-50;
aa(theta9)=1;

hold on; stem(aa)
xlabel('DOA Estimation (degrees)');
ylabel('Spatial Spectrum');

