%{
file name: OC_rock_clay_Li_nondim.m

author: Evan J Ramos
date:   20 Apr 2023

description: this script models nondimensional versions of ODES for OC 
cycling (2 pool) and clay formation/dissolution, rock dissolution, and Li
isotope transfer.

This formulation includes accounting for NPP and its relationship with soil
formation, C accumulation and its relationship with clay content, and clay
formation and its relationship to rock dissolution. Dissolution of clay and
rock is related to OC-enhanced dissolution and scaled to relative
masses and specific surface areas (following literature parameterization of
ligand promoted dissolution).


FULL EQUATION #1 BELOW (slow C degradation FULLY assoc. w minerals)
M(1): dimensionless C_fast
M(2): dimensionless C_slow
M(3): dimensionless X_r
M(4): dimensionless X_cl
M(5): dimensionless Li_cl
M(6): dimensionless Li_w
M(7): d7Li_w

ODEs = @(t,M) [T_Cfast(t)*(1 - M(1)); ... C_fast
               fC*T_Cfast(t)*M(1) - T_Cslow(t)*M(2)*(1/M(4)).^n; ... C_slow 
               -((1-fC)*T_Cfast(t)*m_fast(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + ...   X_r
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3); ...
               fCl*((1-fC)*T_Cfast(t)*m_fast(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + ... X_cl
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3) - ...
                ((1-fC)*T_Cfast(t)*m_fast(M(1),M(2),M(3),M(4))*a_cl(M(3),M(4)) + ...
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(1),M(2),M(3),M(4))*a_cl(M(3),M(4)) + T_Cl(t))*M(4); ...
                fCl*((1-fC)*T_Cfast(t)*m_fast(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + ... Li_cl
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3)/M(4)*(M(6)-M(6));...
                fr*((1-fC)*T_Cfast(t)*m_fast(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + ... Li_w
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3)*(1-fCl*DLi*M(6)) + ...
                fr*((1-fC)*T_Cfast(t)*m_fast(M(1),M(2),M(3),M(4))*a_cl(M(3),M(4)) + ...
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(1),M(2),M(3),M(4))*a_cl(M(3),M(4)) + T_Cl(t))*M(4)*DLi*M(5); ...
                fr/M(6)*... d7Li_w
                (((1-fC)*T_Cfast(t)*m_fast(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + ...
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3)*(d7Li_r-M(7)) + ...
                 DLi*Delta_Li*...
                (((1-fC)*T_Cfast(t)*m_fast(M(1),M(2),M(3),M(4))*a_cl(M(3),M(4)) + ...
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(1),M(2),M(3),M(4))*a_cl(M(3),M(4)) + T_Cl(t))*M(4)*M(5) - ...
                 fCl*((1-fC)*T_Cfast(t)*m_fast(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + ... 
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(1),M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3)*M(6)))];
%}

close all
clear
clc

set(0, 'DefaultAxesFontSize', 15,...
       'DefaultLineLinewidth', 2,...
       'DefaultAxesXColor','k',...
       'DefaultAxesYColor','k')

%% Parameters for OC-promoted dissolution

f_wCf = 0; %fraction of transformed fast-cycling C that enhances weathering
           %1 is all, 0 is none
f_wCs = 0; %fraction of transformed slow-cycling C that enhances weathering
           %1 is all, 0 is none

Beta = 1e-1; %factor which tempers the magnitude of weathering enhancement 
             %for slow and fast C

%% Define characteristic timescales for C and clay cycling

%constants
n        = 1.2; %nonlinear exponent for C degradation
rho_s    = 1.7e6; %(g/m^3) average soil density
rho_r    = 2.7e6; %(g/m^3) average rock density (INITIAL ROCK MASS)
rho_w    = 9.97e5; %(g/m^3) water density

%Constants (that go into dimensionless values)
t0       = 1e4;     %(yr) characteristic time
k_Cl     = 1e-6;    %(yr^-1) clay dissolution rate
P_Cfast  = 1500;     %(g m^-3 yr^-1) net primary productivity
k_Cfast  = 1e-1;    %(yr^-1) C degradation rate (fast pool)
k_Cslow  = 4e-4;    %(yr^-1) C degradation rate (slow pool)
k_R      = 1e-5;    %(yr^-1) rock dissolution rate
fCl      = 0.2;     %(1) fraction of rock dissolution -> clay fm.
fC_trans = 0.5;     %(1) fraction of fast C that is converted to slow C
A_r      = 117;     %(m^2/mol) specific surface area of plagioclase
A_cl     = 4776;    %(m^2/mol) SSA of kaolinite

%Li isotope constants
Li_w     = 0.001; %(ppm) initial fluid Li concentration (> 0 for convergence)
Li_r     = 20;    %(ppm) rock Li concentration
Delta_Li = -17;   %per-mil fractionation factor for clay relative to water
DLi      = 200;   %partition coefficient of Li for clay relative to water
d7Li_r   = 3;     %average rock d7Li value

%%%%%%%%%%%%%%%%%%%%%% CHARACTERISTIC TIMESCALES
%characteristic fast C degradation timescale
T_Cfast = @(t) t0*k_Cfast;

%characteristic slow C degradation timescale
T_Cslow = @(t) t0*k_Cslow;

%characteristic clay dissolution timescale
T_Cl = @(t) t0*k_Cl;

%characteristic rock dissolution timescale
T_R  = @(t) t0*k_R;

%%%%%%%%%%%%%%%%%%%%%% MASS AND AREA SCALING TERMS
R_CX = P_Cfast/(rho_r*k_Cfast); %Ratio of s.s. fast C to initial rock mass
m_fast = @(Cf,Xr,Xcl) Beta*f_wCf*R_CX*Cf./(Xr + Xcl);
m_slow = @(Cs,Xr,Xcl) Beta*f_wCs*R_CX*Cs./(Xr + Xcl);

R_A   = A_r/A_cl; %SSA ratio of rock to clay
a_r   = @(Xr,Xcl) R_A*Xr./(R_A*Xr + Xcl);
a_cl  = @(Xr,Xcl) Xcl./(R_A*Xr + Xcl);

fr    = rho_r/rho_w; %relative mass of rock to water

%% Set initial conditions and time vector

%nondimensional time range and time vector for solve
t_min = 0/t0; t_max = t0/t0;
t_fine = linspace(t_min,t_max,5000);

%initial C, clay
Cfast_i  = 0;
Cslow_i  = 0;
Cl_i     = 1e-3; %(> 0 for convergence)
R_i      = 1-Cl_i;
Li_w_i   = 0.00001;
Li_cl_i  = Li_w/Li_r;
d7Li_w_i = d7Li_r;

%% Solve system of EQs

ODEs = @(t,M) [T_Cfast(t)*(1 - M(1)); ... C_fast
               fC_trans*T_Cfast(t)*M(1) - T_Cslow(t)*M(2)*(1/M(4)).^n; ... C_slow 
               -((1-fC_trans)*T_Cfast(t)*m_fast(M(1),M(3),M(4))*a_r(M(3),M(4)) + ...   X_r
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3); ...
               fCl*((1-fC_trans)*T_Cfast(t)*m_fast(M(1),M(3),M(4))*a_r(M(3),M(4)) + ... X_cl
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3) - ...
                ((1-fC_trans)*T_Cfast(t)*m_fast(M(1),M(3),M(4))*a_cl(M(3),M(4)) + ...
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(2),M(3),M(4))*a_cl(M(3),M(4)) + T_Cl(t))*M(4); ...
                fCl*((1-fC_trans)*T_Cfast(t)*m_fast(M(1),M(3),M(4))*a_r(M(3),M(4)) + ... Li_cl
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3)/M(4)*(M(6)-M(5));...
                fr*((1-fC_trans)*T_Cfast(t)*m_fast(M(1),M(3),M(4))*a_r(M(3),M(4)) + ... Li_w
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3)*(1-fCl*DLi*M(6)) + ...
                fr*((1-fC_trans)*T_Cfast(t)*m_fast(M(1),M(3),M(4))*a_cl(M(3),M(4)) + ...
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(2),M(3),M(4))*a_cl(M(3),M(4)) + T_Cl(t))*M(4)*DLi*M(5); ...
                fr/M(6)*... d7Li_w
                (((1-fC_trans)*T_Cfast(t)*m_fast(M(1),M(3),M(4))*a_r(M(3),M(4)) + ...
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3)*(d7Li_r-M(7)) + ...
                 DLi*Delta_Li*...
                (((1-fC_trans)*T_Cfast(t)*m_fast(M(1),M(3),M(4))*a_cl(M(3),M(4)) + ...
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(2),M(3),M(4))*a_cl(M(3),M(4)) + T_Cl(t))*M(4)*M(5) - ...
                 fCl*((1-fC_trans)*T_Cfast(t)*m_fast(M(1),M(3),M(4))*a_r(M(3),M(4)) + ... 
                 T_Cslow(t)*(1/M(4)).^n*...
                 m_slow(M(2),M(3),M(4))*a_r(M(3),M(4)) + T_R(t))*M(3)*M(6)))];

%solve system of ODEs
soln = ode45(ODEs,[t_min t_max],[Cfast_i Cslow_i R_i Cl_i Li_cl_i Li_w_i d7Li_r]);
Ms = deval(soln,t_fine);

%% Dimensionalized calculations

%Calculate C percent
Cbulk_perc = (Ms(1,:)+Ms(2,:))*P_Cfast/k_Cfast/rho_s*100;
C_fast_perc  = Ms(1,:)*P_Cfast/k_Cfast/rho_s*100;
C_slow_perc = Ms(2,:)*P_Cfast/k_Cfast/rho_s*100;

%Calculate "soil" Li conc. d7Li value
d7Li_clay = Delta_Li + Ms(7,:);
Li_soil   = (Ms(3,:)*rho_r*Li_r + Ms(4,:)*rho_r.*Ms(6,:))./...
            (Ms(3,:)*rho_r + Ms(4,:)*rho_r);
d7Li_soil = (Ms(3,:)*rho_r*Li_r*d7Li_r + Ms(4,:)*rho_r.*Ms(6,:).*d7Li_clay)./...
            (Ms(3,:)*rho_r*Li_r + Ms(4,:)*rho_r.*Ms(6,:));


%% Visualizing model results

%%%%% Non-dimensionalized model predictions over time
figure(1)
subplot(3,2,1) %OC abundances
hold on
plot(t_fine,Ms(1,:),'k-',t_fine,Ms(2,:),'k--');
legend('fast','slow','interpreter','latex','location','northwest')
ylabel('C$^{\prime}$','interpreter','latex')
set(gca,'xticklabel',[],'ticklabelinterpreter','latex')

subplot(3,2,3) %Rock abundance
hold on
plot(t_fine,Ms(3,:),'k-');
ylabel('X$_{r}^{\prime}$','interpreter','latex')
set(gca,'xticklabel',[],'ticklabelinterpreter','latex')

subplot(3,2,5) %Clay abundance
hold on
plot(t_fine,Ms(4,:),'k-')
ylabel('X$_{cl}^{\prime}$','interpreter','latex')
xlabel('t$^{\prime}$','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')

subplot(3,2,2)
hold on
plot(t_fine,Ms(5,:),'k-')
ylabel('Li$_{cl}^{\prime}$','interpreter','latex')
set(gca,'xticklabel',[],'ticklabelinterpreter','latex')

subplot(3,2,4)
hold on
plot(t_fine,Ms(6,:),'k-');
ylabel('Li$_{w}^{\prime}$','interpreter','latex')
set(gca,'xticklabel',[],'ticklabelinterpreter','latex')

subplot(3,2,6)
hold on
plot(t_fine,Ms(7,:)-d7Li_r,'k-')
ylabel('$\delta^7$Li$_{w}-\delta^7$Li$_{r}$','interpreter','latex')
xlabel('t$^{\prime}$','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')

%%%%% Dimensionalized OC vs. d7Li_clay over time
figgy = figure(2);
set(figgy,'Position',[857 272 373 420])
subplot(2,1,1) %OC_slow vs. d7Li
scatter(d7Li_clay-d7Li_r,...
        (Ms(2,:))*P_Cfast/k_Cfast/rho_r*100,10,t_fine*t0/1000,'filled')
ylabel('C$_{\mathrm{slow}}$ (wt \%)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','xticklabel',[])
cc = colorbar('northoutside');
ylabel(cc,'Time (kyr)','Interpreter','latex')
set(cc,'ticklabelinterpreter','latex')
hold on

subplot(2,1,2) %OC slow + OC_fast vs. d7Li
scatter(d7Li_clay-d7Li_r,...
        (Ms(1,:)+Ms(2,:))*P_Cfast/k_Cfast/rho_r*100,10,t_fine*t0/1000,'filled')
xlabel('$\delta^7$Li$_{clay}-\delta^7$Li$_{rock}$','interpreter','latex')
ylabel('C$_{slow} + $C$_{fast}$ (wt \%)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')
hold on

%%%%% Dimensionalized OC vs. d7Li_clay over time (path illustration)
figgy = figure(3);
set(figgy,'Position',[476 463 368 317])
plot(d7Li_clay-d7Li_r,(Ms(2,:))*P_Cfast/k_Cfast/rho_r*100,'k-')
hold on
plot(d7Li_clay(end)-d7Li_r(end),(Ms(2,end))*P_Cfast/k_Cfast/rho_r*100,...
    'ko','markersize',10,'markerfacecolor','w','linewidth',1.25)
plot(d7Li_clay(1)-d7Li_r(1),(Ms(2,1))*P_Cfast/k_Cfast/rho_r*100,...
    's','color',[224 224 224]/255,'markersize',10,...
    'markerfacecolor','k','linewidth',1.25)
xlabel('$\Delta^7$Li$_{clay-rock}$','interpreter','latex')
ylabel('C$_{slow}$ (wt \%)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')

%%%%% d7Li_soil vs. d7Li_clay over time
figgy = figure(4);
set(figgy,'Position',[100 200 368 317])
plot(d7Li_clay,d7Li_soil,'k-')
hold on
plot(d7Li_clay(end),d7Li_soil(end),...
    'ko','markersize',10,'markerfacecolor','w','linewidth',1.25)
plot(d7Li_clay(1),d7Li_soil(1),...
    's','color',[224 224 224]/255,'markersize',10,...
    'markerfacecolor','k','linewidth',1.25)
xlabel('$\delta^7$Li$_{clay}$','interpreter','latex')
ylabel('$\delta^7$Li$_{soil}$','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')
ylim([0 d7Li_r])