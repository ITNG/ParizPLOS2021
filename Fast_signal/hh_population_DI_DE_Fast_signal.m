%% Pariz, etal, PLOS Com. Bio. 2021
% In This code, for each value of inhomogeneity in input current (DI) 
% and delay among populations (DE), the code is simulated once, and 
% then the location of fast pulses is being found. By defining the signal we rerun 
% the code with the same parameters and apply the signal to selected neurons.
% The saved variable will be the spikes in the population with and without 
% signal, RHO0, RHO, respectively, the location of the applied pulse (Pulse_Loc in time step), and the signal
% itself (ISIG). 
% The code is written in MATLAB R2016a by Aref Pariz.

clc
clear
%% Constants
sec=2.5; % in second
dh=0.01; % in ms
range=1000*sec/dh; % TIME STEPS
N=100; % number of neurons per network

% connetions and weights
numnet=2; % Number of networks
inter_c=.1; % Intera layer connectivity probability
neighbor_c=.05; % inter layer connectivity probability
Nep=0.8; % Portion of excitatory neuron 4E:1I
Nip=1-Nep; % number of inhibitory input for each neuron
Ne=Nep*N;   Ni=N-Ne;
con=[1,1;1,1]; % Matrix shows the intra and inter layer connectivity between two populations
v_thr=20; % Threshold for spike
s0=([-59.8977    0.0536    0.5925    0.3192]);

JEE=0.0025;
JEI=0.005;
JIE=0.01;
JII=0.01;
I1=11;I2=11;
I0=[I1*ones(Ne,1);
    (numnet-1)*I2*ones(Ne,1);
    I1*ones(Ni,1);
    (numnet-1)*I2*ones(Ni,1)];
geeIext=0.003; % Synaptic weight dor incomming Poissonian spike train
taur=0.5;taud=3; %% AMPA and GABAz rise and decay time, respectivley
r=2000; % Piosonian spike train rate
Esyn=[0*ones(numnet*Ne,1);-80*ones(numnet*Ni,1)]; % Synaptic reversal potential
ensemble=1;
in_delay=0.5/dh; % Intra delay between all kind of neurons

%% HH model Parameters
C=ones(Ne,1); % Capacitance
gL_bar=0.3;     gK_bar=36;      gNa_bar=120; % mS/cm^2
E_L=-54.387;    E_K=-77;        E_Na=50; % mV

%% MAIN Part
DE=0:0.25:14; % Inter layer delay
DI=-1:0.05:1; % Inhomogeneity on input current of the first population
for out_delay=DE
    for dI0=DI
        display(out_delay)
        display(dI0)
        fname=['data/RES_de',num2str(out_delay),'_dI_',num2str(dI0),'.mat'];
        if ~exist(fname,'file')
            dI=[dI0*ones(Ne,1);
                (numnet-1)*zeros(Ne,1);
                dI0*ones(Ni,1);
                (numnet-1)*zeros(Ni,1)];
            Nsig=[ones(Ne,1);
                (numnet-1)*zeros(Ne,1);
                ones(Ni,1);
                (numnet-1)*zeros(Ni,1)];
            Pulse_Loc=cell(1,ensemble);
            RHO0=cell(1,ensemble);
            RHO=cell(1,ensemble);
            ISIG=cell(1,ensemble);
            
            tic
            for en=1:ensemble
                display(en)
                A=create_A(numnet,N,Nep,inter_c,neighbor_c,con);
                AB=[];ABnum=[];
                for ii=1:numnet*N
                    xnum=numel(find(A(ii,:)));
                    AB(ii,1:xnum)=find(A(ii,:));
                    ABnum(ii)=xnum;
                end
                delay=create_delay(A,numnet,N,Nep,in_delay,out_delay/dh,con);
                J(1:numnet*Ne,1:numnet*Ne)=JEE;
                J(1:numnet*Ne,numnet*Ne+1:numnet*N)=JIE;
                J(numnet*Ne+1:numnet*N,1:numnet*Ne)=JEI;
                J(numnet*Ne+1:numnet*N,numnet*Ne+1:numnet*N)=JII;
                J=J.*A;
                
                vrnd=5*randn(numnet*N,1);
                v=s0(1)*ones(numnet*N,1)+vrnd;
                m=s0(2)*ones(numnet*N,1);
                h=s0(3)*ones(numnet*N,1);
                n=s0(4)*ones(numnet*N,1);
                rho=(zeros(numnet*N,range))==1;
                S1exci=zeros(numnet*N,1);
                S2exci=zeros(numnet*N,1);
                spike_time1=-inf*(ones(numnet*N,1));
                spike_time2=-inf*(ones(numnet*N,1));
                xS=zeros(numnet*N);
                S=zeros(numnet*N,range+1);
                Isyn=(zeros(numnet*N,1));
                Iext=external_current(numnet*N,range,dh,r,taud);
                tic
                for ii=1:range
                    vi=v;
                    ni=n;mi=m;hi=h;
                    v=v+dh*(dI+I0+geeIext*(0-v).*Iext(:,ii)+Isyn-(gK_bar.*(n.^4).*(v-E_K)+gNa_bar.*(m.^3).*h.*(v-E_Na)+gL_bar.*(v-E_L)));
                    m=m+dh*(.1*((v+40)./(1-exp(-.1*(v+40)))).*(1-m)-4*exp(-0.0556*(v+65)).*m);
                    n=n+dh*(0.01*(v+55)./(1-exp(-.1*(v+55))).*(1-n)-0.125*exp(-0.0125*(v+65)).*n);
                    h=h+dh*(0.07*exp(-0.05*(v+65)).*(1-h)-h./(1+exp(-0.1*(v+35))));
                    fired=(vi<v_thr & v>=v_thr); % if v is above threshod and previous v is less than threshod, neuron fired
                    spike_time1(fired)=spike_time2(fired);
                    spike_time2(fired)=ii; % Spike time at synapse
                    S1exci=((taur./taud).^(taur./(taud-taur))-(taur/taud)^(taud/(taud-taur)))\(exp(-(ii*dh-spike_time1*dh)/taur)-exp(-(ii*dh-spike_time1*dh)/taud));
                    S2exci=((taur./taud).^(taur./(taud-taur))-(taur/taud)^(taud/(taud-taur)))\(exp(-(ii*dh-spike_time2*dh)/taur)-exp(-(ii*dh-spike_time2*dh)/taud));
                    S(:,ii)=S1exci+S2exci; % S1+S2;
                    for jj=1:numnet*N % POST
                        for kk=1:ABnum(jj) % PRE
                            if ii>delay(jj,AB(jj,kk))
                                xS(jj,AB(jj,kk))=J(jj,AB(jj,kk))*S(AB(jj,kk),ii-delay(jj,AB(jj,kk)))*(v(jj)-Esyn(AB(jj,kk)));
                            end
                        end
                    end
                    Isyn = sum(xS,2);
                    rho(:,ii)=fired;
                end
                RHO0{en}=rho;
                toc
                %%
                sigma=2;
                fast_rate_conv;
                rr0=rr;
                [peaks,places]=findpeaks(rr0(1,:),'MinPeakProminence',20);
                freq=ceil(mean(sum(rho(1:Ne,:),2)/sec,1));
                seg=50;
                places(1:seg)=[];
                peaks(1:seg)=[];
                ISI=diff(places);
                pulse_amp=0.25;
                pd=fix(2/dh);
                Isig=zeros(1,range);
                phi=linspace(0,2*pi,seg);
                nn=1;n1=0;
                pulse_loc=zeros(1,seg);
                for kk=1:size(places,2)-1
                    n1=n1+1;
                    pulse_loc(n1)=fix(places(kk)+ISI(kk)*phi(nn)/(2*pi));
                    Isig(pulse_loc(n1):pulse_loc(n1)+pd)=pulse_amp;
                    nn=nn+1;
                    if nn==seg
                        nn=1;
                    end
                end
                %%
                v=s0(1)*ones(numnet*N,1)+vrnd;
                m=s0(2)*ones(numnet*N,1);
                h=s0(3)*ones(numnet*N,1);
                n=s0(4)*ones(numnet*N,1);
                rho=(zeros(numnet*N,range))==1;
                S1exci=zeros(numnet*N,1);
                S2exci=zeros(numnet*N,1);
                spike_time1=-inf*(ones(numnet*N,1));
                spike_time2=-inf*(ones(numnet*N,1));
                xS=zeros(numnet*N);
                S=zeros(numnet*N,range+1);
                Isyn=(zeros(numnet*N,1));
                tic
                for ii=1:range
                    vi=v;
                    v=v+dh*(dI+Isig(ii)*Nsig+I0+geeIext*(0-v).*Iext(:,ii)+Isyn-(gK_bar.*(n.^4).*(v-E_K)+gNa_bar.*(m.^3).*h.*(v-E_Na)+gL_bar.*(v-E_L)));
                    m=m+dh*(.1*((v+40)./(1-exp(-.1*(v+40)))).*(1-m)-4*exp(-0.0556*(v+65)).*m);
                    n=n+dh*(0.01*(v+55)./(1-exp(-.1*(v+55))).*(1-n)-0.125*exp(-0.0125*(v+65)).*n);
                    h=h+dh*(0.07*exp(-0.05*(v+65)).*(1-h)-h./(1+exp(-0.1*(v+35))));
                    
                    fired=(vi<v_thr & v>=v_thr); % if v is above threshod and previous v is less than threshod, neuron fired
                    spike_time1(fired)=spike_time2(fired);
                    spike_time2(fired)=ii; % Spike time at synapse
                    S1exci=((taur./taud).^(taur./(taud-taur))-(taur/taud)^(taud/(taud-taur)))\(exp(-(ii*dh-spike_time1*dh)/taur)-exp(-(ii*dh-spike_time1*dh)/taud));
                    S2exci=((taur./taud).^(taur./(taud-taur))-(taur/taud)^(taud/(taud-taur)))\(exp(-(ii*dh-spike_time2*dh)/taur)-exp(-(ii*dh-spike_time2*dh)/taud));
                    S(:,ii)=S1exci+S2exci; % S1+S2;
                    for jj=1:numnet*N % POST
                        for kk=1:ABnum(jj) % PRE
                            if ii>delay(jj,AB(jj,kk))
                                xS(jj,AB(jj,kk))=J(jj,AB(jj,kk))*S(AB(jj,kk),ii-delay(jj,AB(jj,kk)))*(v(jj)-Esyn(AB(jj,kk)));
                            end
                        end
                    end
                    Isyn = sum(xS,2);
                    rho(:,ii)=fired;
                end
                toc
                RHO{en}=rho;
                Pulse_Loc{en}=pulse_loc;
                ISIG{en}=Isig;
            end
            toc
            save(fname,'RHO','RHO0','Pulse_Loc','ISIG');
        end
    end
end
