function   [Pauth_or_an,Pauth_or_sim,Pdet_or_an,Pdet_or_sim,Pauth_and_an,Pauth_and_sim,Pdet_and_an, Pdet_and_sim, ...
    Pauth_csi_an,Pauth_csi_sim, Pdet_csi_an,Pdet_csi_sim, Pauth_cfo_an, Pauth_cfo_sim, Pdet_cfo_an, Pdet_cfo_sim, ...
     Pauth_int_an, Pauth_int_sim, Pdet_int_an,  Pdet_int_sim,mycorr_A,mycorr_E]...
    =csi_cfo_auth_normalized_trad(Ku,thsh_csi,thsh_cfo,f,tau,sampling_rate,Ls,Ns,SNR_AB,SNR_EB,tol_A,tol_B,tol_E,v_AB,intv)

%%%% no combining of metrics here : we use voting, along with the 


CSI_max=1;
CFO_max=1;

thsh_csi=thsh_csi/CSI_max;
thsh_cfo=thsh_cfo/CFO_max;

c=3e8;
rho=besselj(0,2*pi*v_AB*f/c*tau);
noise_init=sqrt(1/(2.*SNR_AB)).*(randn(1,Ku)+1i*randn(1,Ku)); %% noise added to the observation coming from Alice

h_AB_clean=1/sqrt(2).*(randn(1,Ku)+1i*randn(1,Ku)); %% old channel observation from Alice (reference one)

h_AB=h_AB_clean+noise_init; %% old channel observation from Alice (reference one)

noise_A=sqrt(1/(2.*SNR_AB)).*(randn(1,Ku)+1i*randn(1,Ku)); %% noise added to the observation coming from Alice

noise_E=sqrt(1/(2.*SNR_EB)).*(randn(1,Ku)+1i*randn(1,Ku)); %% noise added to the observation coming from Eve

h_AB_new_clean=rho.*h_AB_clean+sqrt(1-rho.^2).*1/sqrt(2).*(randn(1,Ku)+1i*randn(1,Ku));
h_AB_new=h_AB_new_clean+noise_A; %% new observation coming from Alice
h_EB_new_clean=1/sqrt(2).*(randn(1,Ku)+1i*randn(1,Ku)); %% new observation coming from Eve
h_EB_new=h_EB_new_clean+noise_E;

%% CFO estimation

noise_pow_lin_AB=1./SNR_AB; % total noise power
sigma_sq_AB=noise_pow_lin_AB/2;  % half the noise power for the I and Q components of the signal

noise_pow_lin_EB=1./SNR_EB; % total noise power
sigma_sq_EB=noise_pow_lin_EB/2;  % half the noise power for the I and Q components of the signal

Fs=sampling_rate; %% sampling frequency for the CFO estimation

rho_AB=rho;
shift_A=tol_A*f; %% A's  oscillator deviation = (oscillator tolerance * operating frequency)
shift_E=tol_E*f; %% E's  oscillator deviation = (oscillator tolerance * operating frequency)
shift_B=tol_B*f; %% B's  oscillator deviation = (oscillator tolerance * operating frequency)


sig_AB_init=exp(1i*2*pi*(shift_A-shift_B)./Fs.*[0:Ns*Ls-1]').*h_AB_clean; % signal term
sig_AB=exp(1i*2*pi*(shift_A-shift_B)./Fs.*[0:Ns*Ls-1]').*h_AB_new_clean; % signal term 
sig_EB=exp(1i*2*pi*(shift_E-shift_B)./Fs.*[0:Ns*Ls-1]').*h_EB_new_clean; % signal term 
noise_AB_init=sqrt(sigma_sq_AB).*(randn(Ns*Ls,Ku)+1i*randn(Ns*Ls,Ku)); % noise term 
noise_AB=sqrt(sigma_sq_AB).*(randn(Ns*Ls,Ku)+1i*randn(Ns*Ls,Ku)); % noise term 
noise_EB=sqrt(sigma_sq_EB).*(randn(Ns*Ls,Ku)+1i*randn(Ns*Ls,Ku)); % noise term 


y_AB_init=sig_AB_init+noise_AB_init; % received signal 
y_AB=sig_AB+noise_AB; % received signal 
y_EB=sig_EB+noise_EB; % received signal 


% ss=zeros(1,Nsim); %% array containing summations for A-B CFO computation
ss_init=[];
ss=[];
ssE=[];
parfor kk=1:Ku
    %%% A- CFO of the A-B link
  ssp_init=0;
  ssp=0;
  sspE=0;
for ii=1:Ns-1
   ssp_init=ssp_init+(y_AB_init((ii-1)*Ls+1:ii*Ls,kk)'*y_AB_init(ii*Ls+1:(ii+1)*Ls,kk)); 
   ssp=ssp+(y_AB((ii-1)*Ls+1:ii*Ls,kk)'*y_AB(ii*Ls+1:(ii+1)*Ls,kk)); 
   sspE=sspE+(y_EB((ii-1)*Ls+1:ii*Ls,kk)'*y_EB(ii*Ls+1:(ii+1)*Ls,kk)); 
end
ss_init=[ss_init ssp_init];
ss=[ss ssp];
ssE=[ssE sspE];
end

CFO_estim_A_init=1/(2*pi*Ls).*angle(ss_init); % normalized estimated A-B CFO 
CFO_estim_A=1/(2*pi*Ls).*angle(ss); % normalized estimated A-B CFO 
CFO_estim_E=1/(2*pi*Ls).*angle(ssE); % normalized estimated E-B CFO 
TS_cfo_A=abs(CFO_estim_A-CFO_estim_A_init)./CFO_max;
TS_cfo_E=abs(CFO_estim_E-CFO_estim_A_init)./CFO_max;



%%% we simulate here the CFO according to the CRLB approximation (Gaussian
%%% RV where the actual CFO is its mean value)
% CFO_estim_A_init_biss=(shift_A-shift_B)./Fs+1./sqrt(4*pi^2*Ls^3*(Ns-1)*SNR_AB.*abs(h_AB_clean).^2).*randn(1,Ku); % normalized estimated A-B CFO 
% CFO_estim_A_biss=(shift_A-shift_B)./Fs+1./sqrt(4*pi^2*Ls^3*(Ns-1)*SNR_AB.*abs(h_AB_new_clean).^2).*randn(1,Ku); % normalized estimated A-B CFO 
% CFO_estim_E_biss=(shift_E-shift_B)./Fs+1./sqrt(4*pi^2*Ls^3*(Ns-1)*SNR_EB.*abs(h_EB_new_clean).^2).*randn(1,Ku); % normalized estimated E-B CFO 
% TS_cfo_A_biss=abs(CFO_estim_A_biss-CFO_estim_A_init_biss)./CFO_max;
% TS_cfo_E_biss=abs(CFO_estim_E_biss-CFO_estim_A_init_biss)./CFO_max;


%%% test stats for CSI
TS_A=(abs(h_AB_new-h_AB)).^2./CSI_max; %% test statistic of A's observation
TS_E=(abs(h_EB_new-h_AB)).^2./CSI_max; %% test statistic of E's observation




%% 2- Proposed scheme - flipping on 1 interval (working on one of the two intervals)
[ind_A_misc]=find(TS_A>thsh_csi);
[ind_E_misc]=find(TS_E<=thsh_csi);
[ind_A_well]=find(TS_A<=thsh_csi);
[ind_E_well]=find(TS_E>thsh_csi);



[TS_A_misc]=TS_A(find(TS_A>thsh_csi));
[TS_E_misc]=TS_E(find(TS_E<=thsh_csi));
[TS_A_well]=TS_A(find(TS_A<=thsh_csi));
[TS_E_well]=TS_E(find(TS_E>thsh_csi));
L_min=0;
L_max=5;
U_min=0;
U_max=10;
% stp=.1;

% sumss_1=zeros(length(L_min:stp:L_max),length(U_min:stp:U_max)); %% contains nbr of mis-classified
% sumss_2=zeros(length(L_min:stp:L_max),length(U_min:stp:U_max)); %% contains nbr of well-classified
% sumss=zeros(length(L_min:stp:L_max),length(U_min:stp:U_max)); %% contains their difference
% 
% range_L=L_min:stp:L_max;
% range_U=U_min:stp:U_max;
%%% a- selecting suboptimal interval


% for ll=1:length(range_L)
%     for uu=1:length(range_U)
% 
% if (intv==1)
% %%% a- Interval 1
%    sumss(ll,uu)= sum(bitand(TS_E_misc<range_U(uu),TS_E_misc>range_L(ll)))-sum(bitand(TS_E_well<range_U(uu),TS_E_well>range_L(ll)));
% else
% %%% b- Interval 2
%    sumss(ll,uu)= sum(bitand(TS_A_misc<range_U(uu),TS_A_misc>range_L(ll)))-sum(bitand(TS_A_well<range_U(uu),TS_A_well>range_L(ll)));
% end
% 
%     end
% end
% 
% max_diff=max(max(sumss));
% [l_ind,u_ind]=find(sumss==max_diff);
% L_star=L_min+(l_ind(1)-1)*stp
% U_star=U_min+(u_ind(1)-1)*stp

if intv==1
L_star=0
U_star=thsh_csi;
else
L_star=thsh_csi;
U_star=U_max;
end
%%% b- flipping 

% N_well_A_flipped=sum(bitand(TS_A_well>=L_star,TS_A_well<=U_star));
% N_well_E_flipped=sum(bitand(TS_E_well>=L_star,TS_E_well<=U_star));
N_misc_A_flipped=sum(bitand(TS_A_misc>=L_star,TS_A_misc<=U_star));
N_misc_E_flipped=sum(bitand(TS_E_misc>=L_star,TS_E_misc<=U_star));
Pauth_s_2=(sum(bitor(TS_A_well<L_star,TS_A_well>U_star))+N_misc_A_flipped)/Ku;
Pdet_s_2=(sum(bitor(TS_E_well<L_star,TS_E_well>U_star))+N_misc_E_flipped)/Ku;


%% part b - KEEP legitimate ones safe by CFO


TS_A_well_in_interv=[TS_A_well(find(bitand(TS_A_well>L_star,TS_A_well<U_star)))];
TS_E_well_in_interv=[TS_E_well(find(bitand(TS_E_well>L_star,TS_E_well<U_star)))] ;
TS_A_misc_in_interv=[TS_A_misc(find(bitand(TS_A_misc>L_star,TS_A_misc<U_star)))];
TS_E_misc_in_interv=[TS_E_misc(find(bitand(TS_E_misc>L_star,TS_E_misc<U_star)))] ;

set_elem_in_interv_csi=[TS_A(find(bitand(TS_A>L_star,TS_A<U_star))) TS_E(find(bitand(TS_E>L_star,TS_E<U_star)))];
% sz=size(set_elem_in_interv)
ind_A_inside_interv=find(bitand(TS_A>L_star,TS_A<U_star));
ind_E_inside_interv=find(bitand(TS_E>L_star,TS_E<U_star));

set_elem_in_interv=[TS_cfo_A(ind_A_inside_interv) TS_cfo_E(ind_E_inside_interv)];

N=length(set_elem_in_interv);

if (~isempty(set_elem_in_interv))

    
   
selec_elements=[];

for ii=1:N
        
          if(bitor(bitand(set_elem_in_interv(ii)<thsh_cfo,intv==1),bitand(set_elem_in_interv(ii)>thsh_cfo,intv==2))) %% parametrized (both)
%         if(set_elem_in_interv(ii)<thsh_cfo) %% interval 1
 
%         if(set_elem_in_interv(ii)>thsh_cfo) %% interval 2
                selec_elements=[selec_elements set_elem_in_interv_csi(ii)];
          end
end



%%% counting if we made correct selections
% 
% sz1=size(TS_E_well_in_interv)
% sz2=size(selec_elements)
if (~isempty(TS_A_well_in_interv) && ~isempty(selec_elements))
N_correctly_well_sel_A=sum(sum(TS_A_well_in_interv==selec_elements'));
else
    N_correctly_well_sel_A=0;
end

if (~isempty(TS_E_well_in_interv) && ~isempty(selec_elements))
N_correctly_well_sel_E=sum(sum(TS_E_well_in_interv==selec_elements'));
else
    N_correctly_well_sel_E=0;
end

if (~isempty(TS_A_misc_in_interv) && ~isempty(selec_elements))
N_wrongly_misc_sel_A=sum(sum(TS_A_misc_in_interv==selec_elements'));
else
    N_wrongly_misc_sel_A=0;
end


if (~isempty(TS_E_misc_in_interv) && ~isempty(selec_elements))
N_wrongly_misc_sel_E=sum(sum(TS_E_misc_in_interv==selec_elements'));
else
    N_wrongly_misc_sel_E=0;
end

else
    N_correctly_well_sel_A=0;
    N_correctly_well_sel_E=0;
    N_wrongly_misc_sel_A=0;
    N_wrongly_misc_sel_E=0;
end


% N_correctly_well_sel_A
% N_correctly_well_sel_E
% N_wrongly_misc_sel_A
% N_wrongly_misc_sel_E
% 
% N_misc_A_flipped
% N_misc_E_flipped



%%% evaluation


Pauth_int_sim=(sum(bitor(TS_A_well<L_star,TS_A_well>U_star))+N_misc_A_flipped+N_correctly_well_sel_A-N_wrongly_misc_sel_A)/Ku;


Pdet_int_sim=(sum(bitor(TS_E_well<L_star,TS_E_well>U_star))+N_misc_E_flipped+N_correctly_well_sel_E-N_wrongly_misc_sel_E)/Ku;





%% 4 - voting 

%% "OR" rule 

Pauth_or_sim=sum(bitor(TS_A<thsh_csi,TS_cfo_A<thsh_cfo))/Ku;

Pdet_or_sim=sum(bitand(TS_E>=thsh_csi,TS_cfo_E>=thsh_cfo))/Ku;

%% "AND" rule 

Pauth_and_sim=sum(bitand(TS_A<thsh_csi,TS_cfo_A<thsh_cfo))/Ku;

Pdet_and_sim=sum(bitor(TS_E>=thsh_csi,TS_cfo_E>=thsh_cfo))/Ku;


%% 5- unique attr

Pauth_csi_sim=sum(TS_A<thsh_csi)/Ku;
Pdet_csi_sim=sum(TS_E>=thsh_csi)/Ku;


Pauth_cfo_sim=sum(TS_cfo_A<thsh_cfo)/Ku;
Pdet_cfo_sim=sum(TS_cfo_E>=thsh_cfo)/Ku;

mycorr=corrcoef(TS_A,TS_cfo_A);
mycorr_A=mycorr(1,2);
mycorr_E=corrcoef(TS_E,TS_cfo_E);
mycorr_E=mycorr_E(1,2);


% 
% mycorr_biss=corrcoef(TS_A,TS_cfo_A_biss);
% mycorr_biss=mycorr_biss(1,2);
% mycorr_E_biss=corrcoef(TS_E,TS_cfo_E_biss);
% mycorr_E_biss=mycorr_E_biss(1,2);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% B- Analytical Part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%% checking analytical of PDF/CDF of CFO's TS
% TS_cfo_A_biss=TS_cfo_A_biss(find(bitand(abs(h_AB_clean).^2>.5/SNR_AB,abs(h_AB_new_clean).^2>.5/SNR_AB)));
% [nn,xx]=hist(TS_cfo_A_biss,100);

% nn=nn./trapz(xx,nn);
% figure
% plot(xx,nn);

% cc=cumsum(nn);cc=cc./max(cc);

% absic=xx(10)
% prob_sim=nn(10)
% 
% cdf_sim=cc(10)
% kmax=30;
% pmax=kmax;
varr_norm=1/(4*pi^2*Ls^3*(Ns-1));
% pdf_CFO_uncond=0;
% cdf_CFO_uncond=0;
% cs=1/4;
% cv=2/3;
% cw=.62;
% landa=SNR_AB;
% f=@(x,y) 1./landa^2*1/(1-rho_AB^2).*exp(-(x+y)./(landa*(1-rho_AB^2))).*...
% besseli(0,2*rho_AB*sqrt(x.*y)./(landa*(1-rho_AB^2)));
% scaling=integral2(f,1,SNR_AB*10,1,SNR_AB*10);
% epss=1;
% F_triple=@(s,v,w) gammaz(1/2+s).*gammaz(v).*gammaz(w).*gammaz(1-s-v+w).*gammaz(v+w-1).*gammaz(2-s-2*v-w)./...
%     (gammaz(-s).*gammaz(1/2+v).*gammaz(1-v).*gammaz(1/2-v)).*...
%     (thsh_cfo.^2*SNR_AB*(1-rho_AB^2)/(2*CFO_max^2*varr_norm)).^(-s).*rho_AB.^(-2*v).*(1).^(-w);
% 
% pdf_triple=2*(1-rho_AB^2)*sqrt(pi)./thsh_cfo.*1/(2*pi*j)^3*...
%     integral3(F_triple,cs-20*1i,cs+20*1i,cv-20*i,cv+20*i,cw-20*1i,cw+20*1i)
% mysum_1=0;
% mysum_2=0;
% argsss=thsh_cfo^2/(2*varr_norm)
% for p=0:pmax

%%%%%% CDF derived by integrating the PDF derived by integral
% F3=@(s,w) gammaz(-2.*s)./gammaz(1-2.*s).*gammaz(1/2+s)./gammaz(-s).*gammaz(w).*gammaz(1-s+kk+w).*gammaz(w-kk-1).*gammaz(2*(kk+1)-s-w).*...
%     (thsh_cfo.^2*SNR_AB*(1-rho_AB^2)/(2*CFO_max^2*varr_norm)).^(-s).*(1).^(-w);
% 
% F3_biss=@(s,w) gammaz(-2.*s)./gammaz(1-2.*s).*gammaz(1/2+s)./gammaz(-s).*gammaz(w).*gammaz(1-kk-w).*gammaz(kk-s).*gammaz(2*kk+w-s)./...
%     (gammaz(3/2-w-kk).*gammaz(kk+w).*gammaz(kk+w-1/2)).*...
%     (thsh_cfo.^2*SNR_AB*(1-rho_AB^2)/(2*CFO_max^2*varr_norm)).^(-s).*(rho_AB^2).^(w);


%%%%% CDF derived directly

% F3=@(s,w) gammaz(s)./gammaz(3/2-s).*gammaz(w).*gammaz(w-kk-1).*gammaz(3/2+kk-s-w).*gammaz(5/2+2*kk-s-w).*...
%     (thsh_cfo.^2*SNR_AB*(1-rho_AB^2)/(2*varr_norm)).^(-s).*(1).^(-w);
% 
% F3_biss=@(s,w) gammaz(s)./gammaz(3/2-s).*gammaz(w).*gammaz(1-kk-w).*gammaz(1/2+kk-s).*gammaz(1/2+2*kk+w-s)./...
%     (gammaz(3/2-w-kk).*gammaz(kk+w).*gammaz(kk+w-1/2)).*...
%     (thsh_cfo.^2*SNR_AB*(1-rho_AB^2)/(2*varr_norm)).^(-s).*(rho_AB^2).^(w);
% 
% % cdf_CFO_uncond=cdf_CFO_uncond+sqrt(pi*SNR_AB/2)*(1-rho_AB^2)^(3/2)*thsh_cfo/sqrt(varr_norm).*(-1).^kk/factorial(kk)*rho_AB^(2*kk)*...
% %     1/(2*pi*1i)^2*(integral2(F3,cs-80*1i,cs+80*1i,cw-80*1i,cw+80*1i)/(factorial(kk)*gamma(1/2-kk)*gamma(1/2+kk))+...
% %     +integral2(F3_biss,cs-80*1i,cs+80*1i,cw-80*1i,cw+80*1i)/rho_AB^2)




%%%%%%%%%% trunctated CDF integration
% 
% 
% F3=@(v,w) gammaz(v).*gammaz(w).*gammaz(3/2+p-v).*gammaz(v+w-1)./(gammaz(1/2+v).*gammaz(1-v).*gammaz(1/2-v).*(1/2+p).*gammaz(1/2+w+p).*(2*v-p-5/2+w)).*...
%     (rho_AB*epss/(SNR_AB*(1-rho_AB^2))).^(-2*v).*(epss/(SNR_AB*(1-rho_AB^2))).^(-w);
% 
% mysum_1=mysum_1+sqrt(pi/2)*thsh_cfo/((1-rho_AB)^2*SNR_AB^2*sqrt(varr_norm)*(2*pi*1i)^2).*(-1)^p/factorial(p).*(thsh_cfo^2/(2*varr_norm))^p*epss^(p+5/2).*...
%        integral2(F3,cv-20*i,cv+20*i,cw-20*i,cw+20*i);
% 
% 
% for k=0:kmax
% F3_biss=@(v,w) gammaz(v).*gammaz(w).*gammaz(1/2+p+w+k).*gammaz(k-1-v-w)./(gammaz(1/2+v).*gammaz(1-v).*gammaz(1/2-v).*(1/2+p).*gammaz(1/2+w+p).*gammaz(v+w+k).*(3/2+p-v+k)).*...
%     (rho_AB/(SNR_AB*(1-rho_AB^2))).^(-2*v).*(epss/(SNR_AB*(1-rho_AB^2))).^(-w);
% 
%     mysum_2=mysum_2+sqrt(pi/2)*thsh_cfo/((1-rho_AB)^2*SNR_AB^2*sqrt(varr_norm)*(2*pi*1i)^2).*(-1)^p/factorial(p).*(thsh_cfo^2/(2*varr_norm))^p*epss^(-p-1/2).*...
%         (-1)^k/factorial(k)*integral2(F3_biss,cv-20*i,cv+20*i,cw-20*i,cw+20*i);
% 
% end
% 
% cdf_temp=mysum_1-mysum_2
%  
% 
% end
% 
% cdf_CFO_uncond=mysum_1-mysum_2



%%%%%%%%%%%%%%%%%%%%%%%%% Gauss Laguerre Solution %%%%%%%%%%%%%%%%%%%%%%



n=25;

syms t;
T=laguerreL(n,t);

Tp=sym2poly(T);

p=flip(roots(Tp));

ww=p./((n+1)^2*laguerreL(n+1,p).^2);

mean_eve=(shift_E(1)-shift_A(1))./Fs;
%% 1- Alice
G=@(u,v) erf(thsh_cfo*sqrt(SNR_AB*(1-rho_AB^2))/(sqrt(2*varr_norm.*(1./u+1./v)))).*besseli(0,2.*rho_AB.*sqrt(u.*v));

%% 2- Eve
G_E=@(u,v) (erf((thsh_cfo-mean_eve)./sqrt(2*varr_norm.*(1./(SNR_AB*u)+1./(SNR_EB*v))))-...
    erf((-thsh_cfo-mean_eve)./sqrt(2*varr_norm.*(1./(SNR_AB*u)+1./(SNR_EB*v)))));

ssv=0;
sseve=0;

for ii=1:n
    for jj=1:n
   ssv=ssv+ww(ii)*ww(jj)*G(p(ii),p(jj)); 
    sseve=sseve+ww(ii)*ww(jj)*G_E(p(ii),p(jj)); 
    end
end

%% single attribute 

%%% 1 - auth- rate (CSI and CFO)

Pauth_cfo_an=(1-rho_AB^2).*ssv;

Pauth_csi_an=1-igamma(1,thsh_csi./(2.*((1-rho_AB).*1+1/SNR_AB)));

%%% 2- det. rate (CSI and CFO)
% sseve
Pdet_cfo_an=1-1/2*sseve;

Pdet_csi_an=igamma(1,thsh_csi./(2+1/SNR_EB+1/SNR_AB));



%% OR rule (Anal.)

Pauth_or_an=Pauth_cfo_an.*Pauth_csi_an+Pauth_cfo_an.*(1-Pauth_csi_an)+(1-Pauth_cfo_an).*Pauth_csi_an;

Pdet_or_an=Pdet_csi_an*Pdet_cfo_an;



%% AND rule (Anal.)

Pauth_and_an=Pauth_cfo_an.*Pauth_csi_an;

Pdet_and_an=Pdet_cfo_an*Pdet_csi_an+Pdet_cfo_an*(1-Pdet_csi_an)+Pdet_csi_an*(1-Pdet_cfo_an);



%% proposed Scheme 




if intv==1
    
%%% interval 1
Pauth_int_an=Pauth_cfo_an.*Pauth_csi_an;

Pdet_int_an=Pdet_cfo_an*(1-Pdet_csi_an)+Pdet_csi_an;

else
%%% interval 2
% 
Pauth_int_an=Pauth_cfo_an*(1-Pauth_csi_an)+Pauth_csi_an;

Pdet_int_an=Pdet_cfo_an.*Pdet_csi_an;
end






end