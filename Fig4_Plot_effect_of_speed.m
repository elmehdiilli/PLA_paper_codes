clear
clc
close all


f=2.5e9; %% 2.4 GHz (e.g., Wi-Fi)
tau=1e-3;
sampling_rate=2e6;
Ls=32;
Ns=2;
c=3e8;


d_A=150;
d_E=150;
d_0=1;
kappa=2;
G_T=10^(5/10);
G_R=10^(5/10);
PL_A=G_T*G_R.*(c./(4*pi*f*d_0)).^2.*(d_0/d_A).^kappa;
PL_E=G_T*G_R.*(c./(4*pi*f*d_0)).^2.*(d_0/d_E).^kappa;

noise=10^(-20.4)*10*10^6;
power=.1e-3;

SNR_dB=10*log10(power./noise); %% this is power over noise
% SNR_dB=30; %% this is power over noise

% PL_A=10.^([[0:2:30]-SNR_dB]./10);
% PL_E=10.^([[0:2:30]-SNR_dB]./10);

SNR_A=PL_A.*10.^(SNR_dB./10); %% SNR by adding the path-loss term
SNR_E=PL_E.*10.^(SNR_dB./10); %% SNR by adding the path-loss term
10*log10(SNR_A)
10*log10(SNR_E)

Ku=10e3; %% number of Alices and Eves


v_B=0; %% fixed authenticator
v_A=[20:5:70]./3.6; %% v [m/s] = v [km/h] /3.6 
v_E=[20:5:70]./3.6; %% v [m/s] = v [km/h] /3.6 
v_AB=v_A-v_B;



thsh_csi=[1];
thsh_cfo=[2e-4];
% tol_A=1e-6*(rand(1,Ku)+.5);
% tol_E=1e-6*(rand(1,Ku)+2);
tol_A=1e-6*ones(1,Ku);
tol_E=1.5e-6*ones(1,Ku);

tol_B=2.5e-6;
intv=2;

Pauth_or_an=zeros(1,length(v_AB));
Pauth_or_sim=zeros(1,length(v_AB));
Pdet_or_an=zeros(1,length(v_AB));
Pdet_or_sim=zeros(1,length(v_AB));
Pauth_and_an=zeros(1,length(v_AB));
Pauth_and_sim=zeros(1,length(v_AB));
Pdet_and_an=zeros(1,length(v_AB));
Pdet_and_sim=zeros(1,length(v_AB));

Pauth_csi_an=zeros(1,length(v_AB));
Pauth_csi_sim=zeros(1,length(v_AB));
Pdet_csi_an=zeros(1,length(v_AB));
Pdet_csi_sim=zeros(1,length(v_AB));


Pauth_cfo_an=zeros(1,length(v_AB));
Pauth_cfo_sim=zeros(1,length(v_AB));
Pdet_cfo_an=zeros(1,length(v_AB));
Pdet_cfo_sim=zeros(1,length(v_AB));


Pauth_int_an=zeros(1,length(v_AB));
Pauth_int_sim=zeros(1,length(v_AB));
Pdet_int_an=zeros(1,length(v_AB));
Pdet_int_sim=zeros(1,length(v_AB));


for ii=1:length(v_AB)
    
[Pauth_or_an(ii),Pauth_or_sim(ii),Pdet_or_an(ii),Pdet_or_sim(ii),Pauth_and_an(ii),Pauth_and_sim(ii),Pdet_and_an(ii), Pdet_and_sim(ii), ...
    Pauth_csi_an(ii),Pauth_csi_sim(ii), Pdet_csi_an(ii),Pdet_csi_sim(ii), Pauth_cfo_an(ii), Pauth_cfo_sim(ii), Pdet_cfo_an(ii), Pdet_cfo_sim(ii),...
     Pauth_int_an(ii), Pauth_int_sim(ii), Pdet_int_an(ii),  Pdet_int_sim(ii)]=...
      csi_cfo_auth_normalized_trad(Ku,thsh_csi,thsh_cfo,f,tau,sampling_rate,Ls,Ns,SNR_A,SNR_A,tol_A,tol_B,tol_E,v_AB(ii),intv) %% oncly CFO
 
end



%%%%%%%%%%%%%%%%%%% vs SNR $$$$$$$$$$$$$$$$$




    figure
    
    
    
    plot(v_AB*3.6,Pauth_csi_an,'k>-');
    hold on
    plot(v_AB*3.6,Pauth_csi_sim,'rx');


    plot(v_AB*3.6,Pauth_cfo_an,'kd-');
    plot(v_AB*3.6,Pauth_cfo_sim,'rx');

    
    
    plot(v_AB*3.6,Pauth_and_an,'ko-');
    plot(v_AB*3.6,Pauth_and_sim,'rx');
    
    plot(v_AB*3.6,Pauth_or_an,'kh-');
    plot(v_AB*3.6,Pauth_or_sim,'rx');
    
    legend('CSI','','CFO','','Strategy I','','Strategy II','Simulation')


    xlabel('$v_{A_iB}$ [km/h]','Interpreter','latex');
    ylabel('$P_{det}$ (Det. Prob.)','Interpreter','latex')

    
  
