
%In this code we assume that the user add in this path the TFSAP toolboox.
%Please cite the following reference:

clc
clear
close all
SampFreq = 256;
delta=4;
num=2;
%addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;
num=4;
Sig1 = 1*exp(1i*(2*pi*(50*t.^2))+1i*(2*pi*(0*t))); %300t或者150t
Sig2 = 1*exp(1i*(2*pi*(50*t.^2))+1i*(2*pi*(30*t))); %300t或者150t
Sig3 = exp(1i*(2*pi*(105*t -35*t.^3)));
Sig4 =1*exp(1i*(2*pi*(120*t -35*t.^3)));
Sig =1*Sig1 +1*Sig4 +Sig3+0.5*Sig2;
IF_O(:,1)=100*t;
IF_O(:,2)=100*t+30;
IF_O(:,3)=-105*t^2+105;
IF_O(:,4)=-100*t^2+120;

IF_O=IF_O/(SampFreq/2);
num=4;
    Sig=awgn(Sig,10,'measured');
% ORIGINAL
fidexmult= extridge_mult_new_modified1(Sig(1:256), num,5);
figure; plot(t,fidexmult)
orderIF = 50; % use the Fourier model to smooth the detected ridge curves；orderIF1 could be larger than the following orderIF
bw = SampFreq/60;%
%bw=10;
delta=5;
alpha = 5;
[extr_Sig1,fidexmult,tfdv] =extridge_mult_new_spec(Sig, SampFreq, num, delta, orderIF,bw,alpha);
thrf = 15;
tic
[findex,interset] = RPRG(fidexmult,thrf);
toc
figure; plot(t,fidexmult)
