clear
clc
close all



f=2.5e9; %% 2.4 GHz (e.g., Wi-Fi)
tau=1e-3;
sampling_rate=2e6;
Ls=32;
Ns=2;
c=3e8;


d_A=500;
d_E=500;
d_0=1;
kappa=2;
G_T=10^(5/10);
G_R=10^(5/10);
PL_A=G_T*G_R.*(c./(4*pi*f*d_0)).^2.*(d_0/d_A).^kappa;
PL_E=G_T*G_R.*(c./(4*pi*f*d_0)).^2.*(d_0/d_E).^kappa;
noise=10^(-20.4)*10*10^6;
power=.5e-3;

SNR_dB=10*log10(power/noise) %% this is power over noise
% SNR_dB=30; %% this is power over noise
% SNR_A=PL_A.*10.^(SNR_dB./10); %% SNR by adding the path-loss term
% SNR_E=PL_E.*10.^(SNR_dB./10); %% SNR by adding the path-loss term

SNR_A=10.^(20./10); %% SNR manually
SNR_E=10.^(20./10); %% SNR manually


% 10*log10(SNR_A)
Ku=10e3; %% number of Alices and Eves
% thsh_csi=[1];
% thsh_cfo=[2e-4];
thsh_csi=[.001:.001:.009 .01:.01:.09 .1:.1:5];
thsh_cfo=[.1e-5:.1e-5:.9e-5 1e-5:1e-5:9e-5 1e-4:5e-5:8e-4];



v_B=0; %% fixed authenticator
v_A=100./3.6; %% v [m/s] = v [km/h] /3.6 
v_E=100./3.6; %% v [m/s] = v [km/h] /3.6 
v_AB=v_A-v_B;

% tol_A=1e-6*(rand(1,Ku)+.5);
% tol_E=1e-6*(rand(1,Ku)+1);
tol_A=1e-6*ones(1,Ku);
tol_E=1.5e-6*ones(1,Ku);

tol_B=2.5e-6;
intv=1;

Pauth_or_an=zeros(length(thsh_csi),length(thsh_cfo));
Pauth_or_sim=zeros(length(thsh_csi),length(thsh_cfo));
Pdet_or_an=zeros(length(thsh_csi),length(thsh_cfo));
Pdet_or_sim=zeros(length(thsh_csi),length(thsh_cfo));
Pauth_and_an=zeros(length(thsh_csi),length(thsh_cfo));
Pauth_and_sim=zeros(length(thsh_csi),length(thsh_cfo));
Pdet_and_an=zeros(length(thsh_csi),length(thsh_cfo));
Pdet_and_sim=zeros(length(thsh_csi),length(thsh_cfo));

Pauth_csi_an=zeros(1,length(thsh_csi));
Pauth_csi_sim=zeros(1,length(thsh_csi));
Pdet_csi_an=zeros(1,length(thsh_csi));
Pdet_csi_sim=zeros(1,length(thsh_csi));


Pauth_cfo_an=zeros(1,length(thsh_cfo));
Pauth_cfo_sim=zeros(1,length(thsh_cfo));
Pdet_cfo_an=zeros(1,length(thsh_cfo));
Pdet_cfo_sim=zeros(1,length(thsh_cfo));


Pauth_int_an=zeros(length(thsh_csi),length(thsh_cfo));
Pauth_int_sim=zeros(length(thsh_csi),length(thsh_cfo));
Pdet_int_an=zeros(length(thsh_csi),length(thsh_cfo));
Pdet_int_sim=zeros(length(thsh_csi),length(thsh_cfo));


for ii=1:length(thsh_csi)
%     thsh_csi(ii)
    for jj=1:length(thsh_cfo)
% thsh_cfo(jj)



[Pauth_or_an(ii,jj),Pauth_or_sim(ii,jj),Pdet_or_an(ii,jj),Pdet_or_sim(ii,jj),Pauth_and_an(ii,jj),Pauth_and_sim(ii,jj),Pdet_and_an(ii,jj), Pdet_and_sim(ii,jj), ...
    Pauth_csi_an(ii),Pauth_csi_sim(ii), Pdet_csi_an(ii),Pdet_csi_sim(ii), Pauth_cfo_an(jj), Pauth_cfo_sim(jj), Pdet_cfo_an(jj), Pdet_cfo_sim(jj),...
     Pauth_int_an(ii,jj), Pauth_int_sim(ii,jj), Pdet_int_an(ii,jj),  Pdet_int_sim(ii,jj)]=...
      csi_cfo_auth_normalized_trad(Ku,thsh_csi(ii),thsh_cfo(jj),f,tau,sampling_rate,Ls,Ns,SNR_A,SNR_E,tol_A,tol_B,tol_E,v_AB,intv)
 
    end
end



%%%%%%%%%%%%%%%%%%% Figure vs SNR $$$$$$$$$$$$$$$$$
figure
% 
%  plot(1-Pauth_int_an(:,end),Pdet_int_an(:,end),'bh-');
%     hold on
%     plot(1-Pauth_int_sim(:,end),Pdet_int_sim(:,end),'b*');

  
      plot(1-Pauth_csi_an,Pdet_csi_an,'kh-');
    hold on
    plot(1-Pauth_csi_sim,Pdet_csi_sim,'k*');
    
    
        plot(1-Pauth_cfo_an,Pdet_cfo_an,'m>-');
    hold on

    plot(1-Pauth_cfo_sim,Pdet_cfo_sim,'m*');
    
    plot(1-Pauth_and_an(:,end),Pdet_and_an(:,end),'go-');
    hold on
    plot(1-Pauth_and_sim(end,:),Pdet_and_sim(end,:),'g*');
    
      plot(1-Pauth_or_an(1,:),Pdet_or_an(1,:),'rs-');
    hold on
    plot(1-Pauth_or_sim(:,end),Pdet_or_sim(:,end),'r*');
    
  
    
%    legend('Proposed Scheme (Strategy 1)/"AND"" rule','','Proposed Scheme (Strategy 2)/"OR"" rule','',...
%        'CSI','','CFO','Simulation')
legend('CSI','CFO', 'Proposed Scheme (Strategy 1)/"AND"" rule','Proposed Scheme (Strategy 2)/"OR"" rule')
      

    xlabel('$P_{fa}$ (False Alarm Prob.)','Interpreter','latex');
    ylabel('$P_{det}$ (Det. Prob.)','Interpreter','latex')



