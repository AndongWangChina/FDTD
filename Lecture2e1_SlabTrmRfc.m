clc;
close all;
clear;
%%
%Question:1 foot materials=304.8mm;mu=2; epsilon=6;
%To get the reflection and transmission
%

%%
%Initialize Simulation
% Define units
% Define constants
c_0=299792458;%units:m/s

%%
%%Define Simulation Parameters
% ? Device parameters
mu_r=2;
epsilon_r=6;
n_max=sqrt(mu_r*epsilon_r);%refractive index of slab
n_bc=1;%refractive index of boundary
thickness=304.8*10^-3;%unit:m
% ? Frequency range (fmax)
f_max=10^9;
lambda_min=c_0/f_max/n_max;

% ? Grid parameters (NRES, etc.)
N_lambda=20;

%% Compute Grid Resolution
Dlambda=lambda_min/N_lambda;
Dd=thickness/4;
DzPrim=min(Dlambda,Dd);%Minimum requirement about Z
N_slab=ceil(thickness/DzPrim);%Steps in slab
Dz=thickness/N_slab;%Z step

%% Buidl grid
N_total=N_slab+2*10+3;%2*10 for spaces;3 for TF/SF and reflection/transmission;
% INITIALIZE Field on Grid
Ey = zeros(1,N_total);
Hx = zeros(1,N_total);
% INITIALIZE MATERIALS TO FREE SPACE
ER = ones(1,N_total);
UR = ones(1,N_total);
%%Build Device on Grid
Nz1=13;Nz2=83;%position of slab
UR(Nz1:Nz2) = ones(1,(Nz2-Nz1+1))*mu_r;
ER(Nz1:Nz2) = ones(1,(Nz2-Nz1+1))*epsilon_r;
%%Compute the Time Step
Dt=n_bc*Dz/2/c_0;

%% Compute Source Parameters
tau=1/2/f_max;%
t_0=6*tau;

%% Compute Number of Time Steps
t_prop=n_max*N_total*Dz/c_0;
T=12*tau+5*t_prop;
STEPS=ceil(T/Dt);

%% Compute the Source Functions 
t=[0:(STEPS-1)]*Dt;%time axis
deltat=n_bc*Dz/2/c_0+Dt/2;%total delay between E and H
A=-sqrt(1/1);%amplitude of H field
Esrc=exp(-((t-t_0)/tau).^2);%amplitude of H field
Hsrc=A*exp(-((t-t_0+deltat)/tau).^2);%amplitude of H field
nzsrc=2;%position of source injection

%% Initialize the Fourier Transforms
N_freq = 100;%Frequency numbers
FREQ = linspace(0,1*10^9,N_freq);%Frequency 
K = exp(-1i*2*pi*Dt.*FREQ);%Kernerl
REF = zeros(1,N_freq);
TRN = zeros(1,N_freq);
SRC = zeros(1,N_freq);

%% Initialize Detection parameters
EyR = zeros(1,N_freq);
EyT = zeros(1,N_freq);


%% FDTD main loop
figure;%for E/H field

% COMPUTE UPDATE COEFFICIENTS
mEy = (c_0*Dt) ./ ER;
mHx = (c_0*Dt) ./ UR;
% Initialize boundary
H2=0; H1=0;
E2=0; E1=0;

% FDTD main loop

for T = 1 : STEPS

    % Record H-Field at Boundary
    H2=H1; H1=Hx(1);
    
    % Update H from E
    for nz = 1 : (N_total-1)
        Hx(nz) = Hx(nz) + mHx(nz)*(Ey(nz+1) - Ey(nz))/Dz;
    end
    Hx(N_total) = Hx(N_total) + mHx(N_total)*(E2 - Ey(N_total))/Dz;
    
%     Handle H Source
    Hx(1)=Hx(1)-mHx(1)/Dz*Esrc(T);
    
    % Record E-Field Boundary
    E2=E1; E1=Ey(N_total);
    
    % Update E from H
    Ey(1) = Ey(1) + mEy(1)*(Hx(1) - H2)/Dz;
    for nz = 2 : N_total
        Ey(nz) = Ey(nz) + mEy(nz)*(Hx(nz) - Hx(nz-1))/Dz;
    end
    
    % Handle E Source
    Ey(2) = Ey(2)-mEy(2)/Dz*Hsrc(T);
    
    % Inject E Source
%     Ey(nzsrc) = Ey(nzsrc) + Esrc(T);
%         Ey(nzsrc) =  Esrc(T);
    % Hx(nzsrc-1) = Hx(nzsrc-1) + Hsrc(T);
    

    subplot(2,1,1);
    % Visualize
    plot([1:N_total]*Dz,Ey,'b');hold on;
    plot([1:N_total]*Dz,Hx,'r');hold off
    ylim([-2 2]);
    xlim([0 0.4035]);
    title(['FIELD STEP ',num2str(T),' of ',num2str(STEPS)]);


    % Update Fourier Transforms
    for nf = 1 : N_freq
       % Integrate REF(nf), TRN(nf), and SRC(nf)
        EyR(nf) = EyR(nf) + (K(nf)^T)*Ey(1);
        EyT(nf) = EyT(nf) + (K(nf)^T)*Ey(N_total);
        SRC(nf) = SRC(nf) + (K(nf)^T)*Esrc(T);
    end
    % FINISH FOURIER TRANSFORMS
%     EyR = EyR*Dt;
%     EyT = EyT*Dt;
%     SRC = SRC*Dt;

    % COMPUTE REFLECTANCE
    % AND TRANSMITTANCE
    REF = abs(EyR./SRC).^2;
    TRN = abs(EyT./SRC).^2;
    CON = REF + TRN;
    

    subplot(2,1,2);
    plot(REF,'r');
    hold on;
    plot(TRN,'b');
    plot(CON,'k--')
    hold off
    ylim([0 2]);
    pause(0.00010);
    
    if ((T-447)^2<0.01)
       a=1
    elseif ((T-837)^2<0.01)
       a=2
    elseif (T-1383)^2<0.01
       a=3
    end
end


%% Post-Process the Data






