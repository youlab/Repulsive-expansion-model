clear all
close all
global  Nstar  Kphi2 Smesh expphi2 kn

%% Parameter Statement
BeiShu=1;

L=1;  % length of interval
tmax=15; % time integration
Tmesh=100;
Smesh=100;
N=100+1;   % number of grid points; spatial step size h=L/(N-1)

kn0=1;
Nstar=.2; %% could change
Kphi2=.8;  %% could change
expphi2=4;   %change

%% Define Initial Conditions
R0=L/(N-1);
IL=1:1:Smesh;
IL0=(exp(log(R0+1)/(Smesh-1)*(IL-1))-1)';
Nu0=1;

vec=[IL0;R0;Nu0];

%% ODE solver
kn=kn0/BeiShu.^2;
[T,R]=ode45(@gene,linspace(0,tmax,Tmesh),vec);
Nu=R(:,end);
Rad=R(:,end-1);
Locations=R(:,1:Smesh);
TimeMatrix=T*ones(1,Smesh);

figure
for i=1:5:Smesh
plot(R(:,i),T,'LineWidth',1)
hold on
end
plot(Rad,T,'k','LineWidth',2)
xlim([0 Rad(end)+0.1])
ylim([0 tmax+1])
set(gca,'YTick',[0 5 10 15],'FontSize',15)
set(gca,'XTick',[0 0.5 1 1.5],'FontSize',15)
xlabel('Position','FontSize',15)
ylabel('Time','FontSize',15)
axis square


%% plot trajectories of positions within a colony
figure
for i=20:8:Smesh
    plot(R(:,i),T,'LineWidth',1.5)
    hold on
end
plot(Rad,T,'k','LineWidth',2)
xlim([0 Rad(end)+0.1])
ylim([0 tmax+1])
set(gca,'YTick',[0 5 10 15],'FontSize',15)
set(gca,'XTick',[0 0.5 1 1.5],'FontSize',15)
xlabel('Position','FontSize',15)
ylabel('Time','FontSize',15)
% set(gca,'DataAspectRatio',[1 7 1])
axis square
box off

