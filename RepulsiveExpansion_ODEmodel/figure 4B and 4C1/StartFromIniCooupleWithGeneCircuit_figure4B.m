clear all
close all
global Nstar Kphi Kphi2 Smesh gamma beta Tstar Pstar da dT dL sigma0 kT kL m expphi expphi2 kn

%% Parameter Statement

BeiShu=1;
L=1;  % length of interval
tmax=251; % time integration
Tmesh=200;
Smesh=100;
N=400+1;   % number of grid points; spatial step size h=L/(N-1)

%%%%%%%% the parameters %%%%%%%%%%%%%%%%%
kT =8500;
kL = 1900;
kP = 979;
Tstar = 1276;
Pstar = 781;
da=0.3;
dT=0.3;
dL=0.0144;
decayP=10800;
sigma0=1; %% basic devide rate

gamma=kP*Tstar/decayP;
a0=0; %change
Nstar=0.9;
Kphi=0.1;
expphi=2;
Kphi2=0.1;
expphi2=4;
m = 4;
beta0=81;
kn0=25;

%% Define Initial Conditions
R0=L/(N-1);
IL=linspace(0,R0,Smesh)';
IL0=IL.^2/R0;

T0=ones(Smesh,1)*0.1;%Ce0*0.1;
L0=ones(Smesh,1)*0;
Nu0=1;

vec=[IL0;T0;L0;a0;R0;Nu0];

%% ODE solver
beta=beta0/BeiShu.^2;%2*kA*height*Nu0/Astar; %% could change
kn=kn0/BeiShu.^2*(1+a0);

[T,R]=ode45(@gene_20171012,linspace(0,tmax,Tmesh),vec);
Nu=R(:,end);
Rad=R(:,end-1)*sqrt(kn0);
Locations=R(:,1:Smesh)*sqrt(kn0);
T7=R(:,Smesh+1:2*Smesh);
Lys=R(:,2*Smesh+1:3*Smesh);
P=gamma*T7.*Lys;
AHL=R(:,3*Smesh+1);
TimeMatrix=T*ones(1,Smesh);

tfinal=min(find(Nu-0.0001<0,1),Tmesh)+floor(27*2/BeiShu);

%% code below to generate figure 4B
figure('units','normalized','outerposition',[0 0 .3 1])
subplot(2,1,1)
surf(Locations,TimeMatrix,T7,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
colormap jet
set(gca,'FontSize',15)
% title('T7RNAP')
ylim([0 tfinal])
xlim([0 Rad(end)+0.1])
set(gca,'YTick',[0 50 100 150],'FontSize',20)
ylabel('Time','FontSize',20)
set(gca,'XTick',[0 0.5 1 1.5],'FontSize',20)
xlabel('Distance','FontSize',20)
view([0 90]);
colorbar
% set (gca,'Ydir','reverse')
% axis square

subplot(2,1,2)
surf(Locations,TimeMatrix,Lys,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
colormap jet
set(gca,'FontSize',15)
ylim([0 tfinal])
xlim([0 Rad(end)+.1])
xlabel('Distance','FontSize',20)
ylabel('Time','FontSize',20)
set(gca,'YTick',[0 50 100 150],'FontSize',20)
set(gca,'XTick',[0 0.5 1 1.5],'FontSize',20)
% title('Lysozyme')
view([0 90]);
colorbar
% set (gca,'Ydir','reverse')
% axis square

%% code below to generate figure 4C base case
figure('units','normalized','outerposition',[0 0 .2 1])
subplot(2,1,1)
plot(T,AHL,'Color',[255 128 0]/255,'LineWidth',3)
xlim([0 101])
ylim([0 1])
set(gca,'YTick',[0 .5 1],'FontSize',15)
set(gca,'XTick',[0 50 100 150],'FontSize',15)
xlabel('Time','FontSize',15)
ylabel('AHL','FontSize',15)
subplot(2,1,2)
plot(Locations(tfinal,:),Lys(tfinal,:).*[ones(Smesh-1,1);0]','r','LineWidth',3)
xlim([0 3.1])
ylim([0 .12])
set(gca,'XTick',[0 1 2 3 4],'FontSize',15)
set(gca,'YTick',[0 0.04 .08],'FontSize',15)
xlabel('Position','FontSize',15)
ylabel('Lysozyme','FontSize',15)
mid=find(Lys(tfinal,:)==min(Lys(tfinal,:)));
RW=(Rad(end)-Locations(tfinal,mid))*sqrt(kn0);
