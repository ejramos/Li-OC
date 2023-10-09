%{
file name: clay_bedrock_Li_OC_mixing.m

author: Evan J Ramos
date:   23 March 2023

description: this script generates mixing curves in OC vs. d7Li space to
test the hypothesis of bedrock and a "true clay" endmember describing the
observed arrays of Li isotopes and OC in Hawaii, Rio Bermejo, and the
Little Deschutes.

%}

close all
clear
clc

set(0, 'DefaultAxesFontSize', 15,...
       'DefaultLineLinewidth', 2,...
       'DefaultAxesXColor','k',...
       'DefaultAxesYColor','k')

%% Prescribe constants

%bedrock attributes (one value chosen)
OC_rock   = 0; %(wt %) bedrock OC content
Li_rock   = 20; %(ppm) bedrock Li concentration
d7Li_rock = 3; %(per mil LSVEC) bedrock d7Li value

%true clay attributes
%{
Note 1: vectors must be the same length
Note 2: all data in array are synthetic, but outputs from
        OC_rock_clay_Li_nondim.m can be used as inputs to the arrays below
%}
OC_clay   = [0 1 1.2 3 4]; %(wt %) clay OC content array
Li_clay   = [0.04 .2 1 8 11]; %(ppm) clay Li concentration
d7Li_clay = [-17 -15 -13 -11 -4]; %(per mil LSVEC) clay d7Li value

%fraction clay
f_clay = linspace(0,1,100000)';

%preallocating mixing arrays of Li, d7Li, and OC for each clay endmember
Li_samp   = nan(length(f_clay),length(OC_clay));
d7Li_samp = nan(length(f_clay),length(OC_clay));
OC_samp   = nan(length(f_clay),length(OC_clay));

%% Calculate d7Li_samp and OC_samp for each clay endmember

for i = 1:length(OC_clay)
    Li_samp(:,i)   = f_clay*Li_clay(i) + (1-f_clay)*Li_rock;
    d7Li_samp(:,i) = (f_clay*Li_clay(i)*d7Li_clay(i) + ...
                     (1-f_clay)*Li_rock*d7Li_rock)./Li_samp(:,i);
    OC_samp(:,i)   = f_clay*OC_clay(i) + (1-f_clay)*OC_rock;
end
%% Visualize endmembers and mixing arrays

figure('position',[588 360 448 315])
for i = 1:length(OC_clay)
    scatter(d7Li_samp(:,i),OC_samp(:,i),10,f_clay,'filled',...
            'markeredgecolor','none')
    hold on
end
hold on
plot(d7Li_clay,OC_clay,'k^',...
    'markersize',10,'LineWidth',1,'markerfacecolor','w')
plot(d7Li_rock,OC_rock,'ko',...
    'markersize',10,'LineWidth',1,'markerfacecolor','w')

xlabel('$\delta^7$Li','interpreter','latex')
ylabel('OC$_{slow}$ (\%)','interpreter','latex')
cc = colorbar;
ylabel(cc,'Fraction clay','interpreter','latex')
set(cc,'ticklabelinterpreter','latex')
set(gca,'ticklabelinterpreter','latex')