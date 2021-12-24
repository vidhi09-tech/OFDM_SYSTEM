%OFDM SYSTEM--AND -BER OF QPSK -
clc;
clear all;
save all;
close all;
N = 64; % number of bits or symbols

% Transmitter
symbol=10^3;
k=8; %No of Channel coefficient
n=0:k;
a=0.5;  %exponetial decay coefficient
h=exp(-a*n); %channel coefficient
h1=h/(sqrt(sum(h.*h)));
b=rand(1,k+1);
h1i=h1.*cos(b);
h1q=h1.*sin(b);
hh=h1i+j.*h1q;
hf1=fft(hh,64,2); % Channel in frequency domain
hf=kron(ones(symbol,1),hf1); %Converting to matrix form

       
Eb_N0_dB = [0:10];% multiple Eb/N0 values
for ii = 1:length(Eb_N0_dB)
ip=rand(1,N*symbol)>0.5 ;% generating 0,1 with equal probability
ip1=rand(1,N*symbol)>0.5; % generating 0,1 with equal probability
s = 1/sqrt(2)*[(2*ip-1)+j*(2*ip1-1)]; % QPSK modulation 0 -> -1; 1 -> 0

s1=reshape(s,N,symbol);
y8=ifft(s1);
y9=reshape(y8,1,N*symbol);
y=y8.';
y4=[y(:,[49:64]) y];
pow=sqrt(sum(abs(y9).^2)/(N*symbol));

%Convolution
for i=1:symbol
    y5(i,:)=conv(hh,y4(i,:));
end

%receiver
y7=y5;
y1=reshape(y7.',1,symbol*(80+k));
n= 1/sqrt(4)*[randn(1,symbol*(80+k))+j*randn(1,symbol*(80+k))];
y2=y1+10^(-Eb_N0_dB(ii)/20)*n*pow; % Noise addition
y3=reshape(y2,80+k,symbol).';
y3=y3(:,[17:80]);  %yt = yt(:,[17:80]);
z=fft(y3.');
z2=z.';
z3=z2./hf;
z4=z3.';
z1=reshape(z4,1,symbol*N);
z5=reshape(z,1,symbol*N);

%Hard decision
ipHat1 = real(z1)>0;
ipHat2 = imag(z1)>0;
ipHat3 = real(z5)>0;
ipHat4 = imag(z5)>0;

% counting the errors
nErr1(ii) = size(find([ip- ipHat1]),2);
nErr2(ii) = size(find([ip1- ipHat2]),2);
nErr3(ii) = size(find([ip- ipHat3]),2);
nErr4(ii) = size(find([ip1- ipHat4]),2);
% nErr1(ii) = size(find([ip- ipHat1]),2);
% nErr3(ii) = size(find([ip1- ipHat3]),2);
end

simBer = (nErr1+nErr2)/(2.*N*symbol) % simulated ber
simBerfad=(nErr3+nErr4)/(2.*N*symbol) %fading BER
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))) % theoretical ber
theoryBer1 = 0.5*erfc(sqrt(0.5.*(10.^(Eb_N0_dB/10)))); % theoretical ber
figure
semilogy(Eb_N0_dB,theoryBer,'b+-');
hold on
semilogy(Eb_N0_dB,simBer,'r-*');
hold on
semilogy(Eb_N0_dB,simBerfad,'k');
% hold on
% semilogy(Eb_N0_dB,theoryBer1,'k--');
grid on