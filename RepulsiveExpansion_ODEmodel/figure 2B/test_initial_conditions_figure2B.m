clear all
close all
global  Nstar  Kphi2 Smesh expphi2 kn

%% Parameter Statement
BeiShu=1;

L=1;  % length of interval
tmax=21; % time integration
Tmesh=200;
Smesh=100;
N=100+1;   % number of grid points; spatial step size h=L/(N-1)

kn0=1;
Nstar=.2; %% could change
Kphi2=.8;  %% could change
expphi2=4;   %change
kn=kn0/BeiShu.^2;

%% generate figure2B
figure('units','normalized','outerposition',[0 0 .5 .8])
for N=[10 50 100 200 400]
    %% Define Initial Conditions
    R0=L/(N-1);
    IL=1:1:Smesh;
    IL0=(exp(log(R0+1)/(Smesh-1)*(IL-1))-1)';

    Nu0=1;
    vec=[IL0;R0;Nu0];
 
    %% ODE solver
    [T,R]=ode45(@gene,linspace(0,tmax,Tmesh),vec);
    Nu=R(:,end);
    Rad=R(:,end-1);
    Locations=R(:,1:Smesh);
    TimeMatrix=T*ones(1,Smesh);
    subplot(4,3,1:3)
    plot([IL0; IL0(end)],[ones(Smesh,1)*1; 0],'LineWidth',3)
    hold on
    xlim([0 0.12])
    ylim([0 1.05])
    set(gca,'XTick',[0 0.1 .2 .3 .4 .5],'FontSize',20)
    set(gca,'YTick',[0 1],'FontSize',20)
%     xlabel('Position','FontSize',15)
%     ylabel('Cell Density','FontSize',15)
    subplot(4,3,4:12)
    plot(Rad,T, 'LineWidth',3)
    hold on
    ylim([0 tmax+1])
    xlim([0 1.6])
    set(gca,'YTick',[0 10 20],'FontSize',20)
    set(gca,'XTick',[0 0.5 1 1.5 2],'FontSize',20)
    ylabel('Time','FontSize',20)
    xlabel('Colony Radius','FontSize',15)  
%     set(gca,'DataAspectRatio',[1 7 1])
    box off
end