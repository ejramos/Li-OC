%%
%{
filename: OC_Li_regression.m

author: Evan J. Ramos
date:   05 January 2023

Description: This script performs Deming regressions of OC abundances and 
Li isotope ratios from Hawaii (soil averages), the Little Deschutes river 
in Oregon, USA (river suspended sediment) and the Rio Bermejo, Argentina
(river suspended sediment). 

%}

close all
clear
clc

set(0, 'DefaultAxesFontSize', 15,...
       'DefaultLineLinewidth', 2,...
       'DefaultAxesXColor','k',...
       'DefaultAxesYColor','k')

%% Assign arrays


%{
See excel spreadsheets in this repository for the complete dataset
%}

%%%%%%%%%%%%%%%%%%%%%%% HAWAII
%column 1: Mean annual precipitation (mm)
%column 2: average d7Li value (per mil LSVEC)
%column 3: standard error d7Li value
%column 4: average OC content (wt %)
%column 5: standard error OC content

HAWAII = [2500	3.3	 0.1 6.0  2.6;
          2500	3.6	 0.1 5.3  2.6;
          2500	-0.1 1.0 13.5 2.1;
          2500  2.4  1.5 9.5  2.9;
          2500  2.9  1.3 5.5  5.5;
          2500  4.2  0.8 3.9  2.7;
          180   5.3  0.5 4.1  0.7;
          570   6.7  1.0 7.6  0.9;
          1260  8.2  1.8 9.2  2.9;
          1380  6.4  0.6 10.2 2.9;
          2500  0.7  0.4 15.4 2.4;
          3000  3.2  0.4 26.3 6.2;
          1730  5.6  0.9 2.2  0.8;
          350   10.1 1.6 1.1  0.2];

d7Li_rock_H = 4;

%%%%%%%%%%%%%%%%%%%%%%% LITTLE DESCHUTES
%column 1: d7Li value (per mil LSVEC)
%column 2: OC content (wt %)

LD = [-5.4	8.1;
      -5.4	7.9;
      -6.6	10.6
      -10.1	15.6
      -8.8	14.3
      -7.3	12.8
      -6.5	8.5
      -7.6	16.7
      -6.9	13.8
      -6.4	13.3
      -6.8	9
      -7.6	13.3
      -6.0	10.1
      -5.2	14.6];

d7Li_rock_LD = 3.3;

%%%%%%%%%%%%%%%%%%%%%%% RIO BERMEJO
%column 1: d7Li value (per mil LSVEC)
%column 2: OC content (wt %)
RB = [-0.9	0.31
      -1.2	0.23
      -0.7	0.28
      0.3	0.20
      -0.4	0.22
      -1.9	0.30
      4.1	0.06
      -2.6	0.36
      -1.7	0.27
      -1.6	0.30
      -2.3	0.30
      -1.9	0.27
      -0.9	0.25];

d7Li_rock_RB = 4.1;

%%%%%% error for measurements of Li isotope ratios and OC abundances
Li_err = 0.5; %long-term precision of reference material
OC_err = 0.1; %long-term precision of reference material

%% Perform Deming regressions and plot arrays

%%%%%%%%%%%%%%%%%%%%%%%%%

% Superimposing all arrays


%%%%%%%%%%%%%%%%%% HAWAII ARRAY
%select samples that meet criteria for weathering-OC inquiries (see
%Supporting Information for the justification)
sel = HAWAII(:,1) >= 1700 & HAWAII(:,1) < 3000; %isolating MAP window

figure('Position',[305 325 521 372])
%perform Deming regression
[ p, saxipbi, xplot, low_alpha, high_alpha] = ...
                linfitxy(HAWAII(sel,2)-d7Li_rock_H, HAWAII(sel,4), ...
                         HAWAII(sel,3), HAWAII(sel,5));
x_corr = xplot >= min(HAWAII(sel,2)-d7Li_rock_H) & ...
         xplot <= max(HAWAII(sel,2)-d7Li_rock_H);
xplot = xplot(x_corr);
saxipbi = saxipbi(:,x_corr);

%quantify the p-value for the regression (without error)
[~,~,~,~,stats] = regress(HAWAII(sel,4),...
    [ones(length(HAWAII(sel,4)),1) HAWAII(sel,2)-d7Li_rock_H]);

%print statistics of the regression to command window
fprintf('Hawaii rho = %.3f\n',corr(HAWAII(sel,2)-d7Li_rock_H,HAWAII(sel,4)))
fprintf('Hawaii p-value = %.4f\n\n',stats(3))

%generate the plot for Hawaii array
patch([xplot fliplr(xplot)],...
      [saxipbi(low_alpha,:) fliplr(saxipbi(high_alpha,:))],...
      'k','facealpha',0.2,'linestyle','none')
hold on
plot([min(xplot),max(xplot)], p(1)*[min(xplot),max(xplot)]+p(2),...
     'Color','k','handlevisibility','off')
errorbar(HAWAII(sel,2)-d7Li_rock_H,HAWAII(sel,4),...
         HAWAII(sel,5),HAWAII(sel,5),HAWAII(sel,3),HAWAII(sel,3),...
        'ko','markerfacecolor','w','markersize',10,'capsize',0,...
        'handlevisibility','off')

%%%%%%%%%%%%%%%%%% LITTLE DESCHUTES ARRAY
%perform Deming regression
[ p, saxipbi, xplot, low_alpha, high_alpha] = ...
                linfitxy(LD(:,1)-d7Li_rock_LD, LD(:,2), ...
                         Li_err*ones(length(LD),1), ...
                         OC_err*ones(length(LD),1));
x_corr = xplot >= min(LD(:,1)-d7Li_rock_LD) & ...
         xplot <= max(LD(:,1)-d7Li_rock_LD);
xplot = xplot(x_corr);
saxipbi = saxipbi(:,x_corr);

%quantify p-value for regression (without error)
[~,~,~,~,stats] = regress(LD(:,2),...
                          [ones(length(LD),1) LD(:,1)-d7Li_rock_LD]);

fprintf('Little Deschutes rho = %.3f\n',corr(LD(:,1)-d7Li_rock_LD,LD(:,2)))
fprintf('Little Deschutes p-value = %.4f\n\n',stats(3))

%generate plot for Little Deschutes array
patch([xplot fliplr(xplot)],...
      [saxipbi(low_alpha,:) fliplr(saxipbi(high_alpha,:))],...
      [0 128 128]/255,'facealpha',0.2,'linestyle','none')

plot([min(xplot),max(xplot)], p(1)*[min(xplot),max(xplot)]+p(2),...
     'Color',[0 128 128]/255,'handlevisibility','off')
errorbar(LD(:,1)-d7Li_rock_LD,LD(:,2),...
         OC_err*ones(length(LD),1),OC_err*ones(length(LD),1),...
         Li_err*ones(length(LD),1),Li_err*ones(length(LD),1),...
        '^','markerfacecolor','w','markersize',10,'capsize',0,...
        'handlevisibility','off','color',[0 128 128]/255)

%%%%%%%%%%%%%%%%%% RIO BERMEJO ARRAY
%perform Deming regression
[ p, saxipbi, xplot, low_alpha, high_alpha] = ...
                linfitxy(RB(:,1)-d7Li_rock_RB, RB(:,2), ...
                         Li_err*ones(length(RB),1), ...
                         OC_err*ones(length(RB),1));
x_corr = xplot >= min(RB(:,1)-d7Li_rock_RB) & ...
         xplot <= max(RB(:,1)-d7Li_rock_RB);
xplot = xplot(x_corr);
saxipbi = saxipbi(:,x_corr);

%quantify p-value for regression (without error)
[~,~,~,~,stats] = regress(RB(:,2),...
                          [ones(length(RB),1) RB(:,1)-4.1]);

fprintf('Rio Bermejo rho = %.3f\n',corr(RB(:,1)-d7Li_rock_RB,RB(:,2)))
fprintf('Rio Bermejo p-value = %.4f\n\n',stats(3))

%generate plot for Rio Bermejo array
patch([xplot fliplr(xplot)],...
      [saxipbi(low_alpha,:) fliplr(saxipbi(high_alpha,:))],...
      [128 128 0]/255,'facealpha',0.2,'linestyle','none')

plot([min(xplot),max(xplot)], p(1)*[min(xplot),max(xplot)]+p(2),...
     'Color',[128 128 0]/255,'handlevisibility','off')
errorbar(RB(:,1)-d7Li_rock_RB,RB(:,2),...
         OC_err*ones(length(RB),1),OC_err*ones(length(RB),1),...
         Li_err*ones(length(RB),1),Li_err*ones(length(RB),1),...
        'v','markerfacecolor','w','markersize',10,'capsize',0,...
        'handlevisibility','off','color',[128 128 0]/255)
xlim([-16 2])
ylim([0 20])
set(gca,'ticklabelinterpreter','latex')
xlabel('$\Delta^7$Li$_{fine-rock}$','interpreter','latex')
ylabel('OC (\%)','interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot for solely the Rio Bermejo array


figure('Position',[888 91 499 356])
patch([xplot fliplr(xplot)],...
      [saxipbi(low_alpha,:) fliplr(saxipbi(high_alpha,:))],...
      [128 128 0]/255,'facealpha',0.2,'linestyle','none')
hold on
plot([min(xplot),max(xplot)], p(1)*[min(xplot),max(xplot)]+p(2),...
     'Color',[128 128 0]/255,'handlevisibility','off')

errorbar(RB(:,1)-d7Li_rock_RB,RB(:,2),...
         OC_err*ones(length(RB),1),OC_err*ones(length(RB),1),...
         Li_err*ones(length(RB),1),Li_err*ones(length(RB),1),...
        'v','markerfacecolor','w','markersize',10,'capsize',0,...
        'handlevisibility','off','color',[128 128 0]/255)
xlim([-8 0])
ylim([0 .5])
xlabel('$\Delta^7$Li$_{fine-rock}$','interpreter','latex')
ylabel('OC (\%)','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')
