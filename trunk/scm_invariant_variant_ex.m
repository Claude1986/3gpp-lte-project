% SCM 
%_SCM 3GPP Spatial Channel Model 
% % SCM-2006-07-05 V1.21
% MATLAB implementation
% of the 3GPP Spatial Channel Model
% (3GPP TR 25.996)
% Implementation Documentation
% Version: 1.21 (scm-05-07-2006)
% Date: July 5, 2006
% File: scm-05-07-2006.zip
clear all
close all
clc;
Nfft=64;  % FFT size
L=6;   % # of taps Tap scmparset içinde
CP=Nfft/4;
F=dftmtx(Nfft);
FL=F(:,1:L); 
t=Nfft+CP;  %N+CP ==> Length of the OFDM symbol 
V = [0.1 3 60 120]  %[km/h] kullanýcýlarýn hýzý

for di=1:3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Channel Parameter Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath winner  % adds directory for the MIMO channel 
% initialization
scmpar=scmparset;
linkpar=linkparset(1); % linkpar=linkparset(1) <==> linkpar=linkparset
antpar=antparset;
%Configuring parameter to generate impulse response for each OFDM symbol. 
SlotDuration=0.5e-3; %0.5 ms
N.OFDM_symbol_pr_slot=6; % 6 symbols pr. slot
speed_of_light=2.99792458e8; %m/s
fcarrier=2e9; % 2GHz carrier frequency
wavelength=speed_of_light/fcarrier;
delta_t=SlotDuration/N.OFDM_symbol_pr_slot; %OFDM symbol duration
linkpar.MsVelocity=1000/3600*V(di); % unit [m/sec] Velocity of UE is configured
%Configuring parameter to generate impulse response for each OFDM symbol.
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
scmpar.SampleDensity = wavelength / (linkpar.MsVelocity*2*delta_t);
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
scmpar.NumTimeSamples=t; % generata scm for each ofdm symbol duration  
scmpar.XpdIndependentPower='yes'; % Power Delay Profile normalization
scmpar.ScmOptions='los'; % Power Delay Profile normalization

[H1 delays out]=scm(scmpar,linkpar,antpar);
%%%%!!!!!!!!!!!!!! end of channel matrix generation%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Channel matrix generation  using final conditions as initial conditions in next function call
%for i=1:2
[H1 delays out]=scm(scmpar,linkpar,antpar);
[H2 delays out]=scm(scmpar,linkpar,antpar,out);

h11=squeeze(H1(1,1,:,:)); % channel taps
H11=FL*h11; % Channel frequency. response

h11_2=squeeze(H2(1,1,:,:));
H11_2=FL*h11_2;

h12=squeeze(H1(1,2,:,:));
H12=FL*h12;
h21=squeeze(H1(2,1,:,:));
H21=FL*h21;
h22=squeeze(H2(2,2,:,:));
H22=FL*h22;


figure
subplot(2,1,1)
surf(abs(H11))
title(['Kanal Transfer Fonksiyonu |H_{11}|, V= '  num2str(V(di))  ' km/h'])
xlabel('Zaman [n]')
ylabel('Frekans [k]')
zlabel('|H_{11}[k,n]|')
%axis([0 t 0 N zmin zmax])

% subplot(2,2,2)
% surf(abs([H11 H11_2]))
% title(['H_{11}, cont. channel [H11(1) H11(2)], V= '  num2str(V(di))  ' km/h'])
% xlabel('time')
% ylabel('Frekans')
% zlabel('Mutlak büyüklük')

subplot(2,1,2)
contourf(abs(H11))
colorbar
title(['Kanal Transfer Fonksiyonu |H_{11}|, V= '  num2str(V(di))  ' km/h'])
xlabel('Zaman [n]')
ylabel('Frekans [k]')
zlabel('|H_{11}[k,n]|')
%axis([0 t 0 N zmin zmax])
% 
% subplot(2,2,4)
% contourf(abs([H11 H11_2]))
% colorbar
% title(['H_{11}, cont. channel [H11(1) H11(2)], V= '  num2str(V(di))  ' km/h'])
% xlabel('time')
% ylabel('Frekans')
% zlabel('Mutlak büyüklük')
% saveas(gcf,['H_{11}' num2str(di)],'fig')

figure

subplot(2,2,1)
contourf(abs(H11))
colorbar
title(['Kanal Transfer Fonksiyonu |H_{11}|, V= '  num2str(V(di))  ' km/h'])
xlabel('Zaman [n]')
ylabel('Frekans [k]')
zlabel('|H_{11}[k,n]|')

subplot(2,2,2)
contourf(abs(H12))
colorbar
title(['Kanal Transfer Fonksiyonu |H_{12}|, V= '  num2str(V(di))  ' km/h'])
xlabel('Zaman [n]')
ylabel('Frekans [k]')
zlabel('|H_{12}[k,n]|')

subplot(2,2,3)
contourf(abs(H21))
colorbar
title(['Kanal Transfer Fonksiyonu |H_{21}|, V= '  num2str(V(di))  ' km/h'])
xlabel('Zaman [n]')
ylabel('Frekans [k]')
zlabel('|H_{21}[k,n]|')
subplot(2,2,4)
contourf(abs(H22))
colorbar
title(['Kanal Transfer Fonksiyonu |H_{22}|, V= '  num2str(V(di))  ' km/h'])
xlabel('Zaman [n]')
ylabel('Frekans [k]')
zlabel('|H_{22}[k,n]|')
saveas(gcf,['H_all' num2str(di)],'fig')


figure
subplot(2,1,1)
contourf(abs(h11))
colorbar
title(['h_{11}(), V= '  num2str(V(di))  ' km/h'])
xlabel('time')
ylabel('taps')
zlabel('Mutlak büyüklük')
%axis([0 t 0 N zmin zmax])

subplot(2,1,2)
contourf(abs(h22))
colorbar
title(['h_{22}, V= '  num2str(V(di))  ' km/h'])
xlabel('time')
ylabel('taps')
zlabel('Mutlak büyüklük')

saveas(gcf,'tapgains','fig')
hvar=convmtx_variant(h11,Nfft+CP);

figure
contourf(abs(hvar))
colorbar
title(['h_{11} conv mtx, V= '  num2str(V(di))  ' km/h'])
xlabel('time index')
ylabel('convolution index')
zlabel('Mutlak büyüklük')
saveas(gcf,['h_{11}convolutionmtx' num2str(di)],'fig')
end
%%


 %% Show that channel tap power gains are unity 
 scmpar.NumTimeSamples=1e5;
[H4 delays out]=scm(scmpar,linkpar,antpar);
mean(abs(sum(H4,3)).^2,4)

%% Autocorrelation of the taps
t_samp=100
N_corr=200


scmpar=scmparset;
linkpar=linkparset(1); % linkpar=linkparset(1) <==> linkpar=linkparset
antpar=antparset;
%Configuring parameter to generate impulse response for each OFDM symbol. 
SlotDuration=0.5e-3; %0.5 ms
N.OFDM_symbol_pr_slot=6; % 6 symbols pr. slot
speed_of_light=2.99792458e8; %m/s
fcarrier=2e9; % 2GHz carrier frequency
wavelength=speed_of_light/fcarrier;
delta_t=SlotDuration/N.OFDM_symbol_pr_slot; %OFDM symbol duration
linkpar.MsVelocity=1000/3600*V(1); % unit [m/sec] Velocity of UE is configured
%Configuring parameter to generate impulse response for each OFDM symbol.
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
scmpar.SampleDensity = wavelength / (linkpar.MsVelocity*2*delta_t);
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
scmpar.NumTimeSamples=t; % generata scm for each ofdm symbol duration  
scmpar.XpdIndependentPower='yes'; % Power Delay Profile normalization
scmpar.ScmOptions='los'; % Power Delay Profile normalization



scmpar.NumTimeSamples=t_samp;
h_tap=zeros(t_samp,N_corr);
corr_mtx=zeros(2*t_samp-1,N_corr);










for i=1:N_corr
[H4 delays out]=scm(scmpar,linkpar,antpar);   
h_tap(:,i)=squeeze(H4(1,1,1,:));
corr_mtx(:,i)=xcorr(h_tap(:,i),'coeff');
end

[corr_res,LAGS]=xcorr(h_tap,'coeff');
figure
mean_corr_res=mean(corr_res,2);
mean_corr_res=mean_corr_res/max(mean_corr_res);
plot(LAGS,abs(mean_corr_res))
grid on
xlabel('Gecikme[n]')
ylabel('R[n]')
saveas(gcf,['Channel autocorrelation mtx'],'fig')
%% Channel Correlation Matrices
N_samp=100
N_frame=2000
scmpar.NumTimeSamples=N_samp;
hTX_taps1=zeros(N_frame*N_samp,2);
hTX_taps2=zeros(N_frame*N_samp,2);
hTX_taps3=zeros(N_frame*N_samp,2);
hTX_taps4=zeros(N_frame*N_samp,2);
hTX_taps5=zeros(N_frame*N_samp,2);
hTX_taps6=zeros(N_frame*N_samp,2);
%--------------RX---
hRX_taps1=zeros(N_frame*N_samp,2);
hRX_taps2=zeros(N_frame*N_samp,2);
hRX_taps3=zeros(N_frame*N_samp,2);
hRX_taps4=zeros(N_frame*N_samp,2);
hRX_taps5=zeros(N_frame*N_samp,2);
hRX_taps6=zeros(N_frame*N_samp,2);






for iFrames=1:N_frame
    [H4 delays out]=scm(scmpar,linkpar,antpar);
    idx = (1:N_samp)+(iFrames-1)*N_samp;
    for it=1:2
        hTX_taps1(idx,it)=squeeze(H4(1,it,1,:));
        hTX_taps2(idx,it)=squeeze(H4(1,it,2,:));
        hTX_taps3(idx,it)=squeeze(H4(1,it,3,:));
        hTX_taps4(idx,it)=squeeze(H4(1,it,4,:));
        hTX_taps5(idx,it)=squeeze(H4(1,it,5,:));
        hTX_taps6(idx,it)=squeeze(H4(1,it,6,:));
        %-----------------RX------------------
        hRX_taps1(idx,it)=squeeze(H4(it,1,1,:));
        hRX_taps2(idx,it)=squeeze(H4(it,1,2,:));
        hRX_taps3(idx,it)=squeeze(H4(it,1,3,:));
        hRX_taps4(idx,it)=squeeze(H4(it,1,4,:));
        hRX_taps5(idx,it)=squeeze(H4(it,1,5,:));
        hRX_taps6(idx,it)=squeeze(H4(it,1,6,:));
    end
end

TxCorrMatrixPath1=corrcoef(hTX_taps1(:,1),hTX_taps1(:,2)).';
TxCorrMatrixPath2=corrcoef(hTX_taps2(:,1),hTX_taps2(:,2)).';
TxCorrMatrixPath3=corrcoef(hTX_taps3(:,1),hTX_taps3(:,2)).';
TxCorrMatrixPath4=corrcoef(hTX_taps4(:,1),hTX_taps4(:,2)).';
TxCorrMatrixPath5=corrcoef(hTX_taps5(:,1),hTX_taps5(:,2)).';
TxCorrMatrixPath6=corrcoef(hTX_taps6(:,1),hTX_taps6(:,2)).';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
RxCorrMatrixPath1=corrcoef(hRX_taps1(:,1),hRX_taps1(:,2)).';
RxCorrMatrixPath2=corrcoef(hRX_taps2(:,1),hRX_taps2(:,2)).';
RxCorrMatrixPath3=corrcoef(hRX_taps3(:,1),hRX_taps3(:,2)).';
RxCorrMatrixPath4=corrcoef(hRX_taps4(:,1),hRX_taps4(:,2)).';
RxCorrMatrixPath5=corrcoef(hRX_taps5(:,1),hRX_taps5(:,2)).';
RxCorrMatrixPath6=corrcoef(hRX_taps6(:,1),hRX_taps6(:,2)).';

TxCorrMatrixPath=zeros(2,2,6);
TxCorrMatrixPath(:,:,1)=corrcoef(hTX_taps1(:,1),hTX_taps1(:,2)).';
TxCorrMatrixPath(:,:,2)=corrcoef(hTX_taps2(:,1),hTX_taps2(:,2)).';
TxCorrMatrixPath(:,:,3)=corrcoef(hTX_taps3(:,1),hTX_taps3(:,2)).';
TxCorrMatrixPath(:,:,4)=corrcoef(hTX_taps4(:,1),hTX_taps4(:,2)).';
TxCorrMatrixPath(:,:,5)=corrcoef(hTX_taps5(:,1),hTX_taps5(:,2)).';
TxCorrMatrixPath(:,:,6)=corrcoef(hTX_taps6(:,1),hTX_taps6(:,2)).';
TxCorr_mean=sum(TxCorrMatrixPath,3)/6

RxCorrMatrixPath=zeros(2,2,6);
RxCorrMatrixPath(:,:,1)=corrcoef(hRX_taps1(:,1),hRX_taps1(:,2)).';
RxCorrMatrixPath(:,:,2)=corrcoef(hRX_taps2(:,1),hRX_taps2(:,2)).';
RxCorrMatrixPath(:,:,3)=corrcoef(hRX_taps3(:,1),hRX_taps3(:,2)).';
RxCorrMatrixPath(:,:,4)=corrcoef(hRX_taps4(:,1),hRX_taps4(:,2)).';
RxCorrMatrixPath(:,:,5)=corrcoef(hRX_taps5(:,1),hRX_taps5(:,2)).';
RxCorrMatrixPath(:,:,6)=corrcoef(hRX_taps6(:,1),hRX_taps6(:,2)).';
RxCorr_mean=sum(RxCorrMatrixPath,3)/6