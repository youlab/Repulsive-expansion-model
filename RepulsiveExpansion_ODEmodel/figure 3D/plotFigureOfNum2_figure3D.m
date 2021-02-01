clear all
close all
global Nstar Kphi Kphi2 Smesh gamma beta Tstar Pstar da dT dL sigma0 kT kL m expphi expphi2 kn

% load SCORE_uniformCell.mat 
BeiShu=1;
Num=108; % 917 877 869 797 660 591 583 584 589 433 466 330 352 202 216 245 108 129-2r
%% Parameter Statement
L=1;  % length of interval
tmax=151; % time integration
Tmesh=200;
Smesh=100;
N=400+1;   % number of grid points; spatial step size h=L/(N-1)

kT =3907.9;
kL =1411.2;
kP = 267.53;
Tstar = 443.21;
Pstar = 435.10;
da=0.3;
dT=0.3;
dL=0.0144;
decayP=10800;
sigma0=1; %% basic devide rate

gamma=kP*Tstar/decayP;
a0 = 0.6;   %change
Nstar= 0.4249;
Kphi = 0.1192;
expphi=6;
Kphi2 =0.7064;
expphi2=3;
m = 4;
beta0=25.728;
kn0  =14.878;

%% Define Initial Conditions
R0=L/(N-1);
IL=1:1:Smesh;
IL0=(exp(log(R0+1)/(Smesh-1)*(IL-1))-1)';

T0=ones(Smesh,1)*0.1;  %Ce0*0.1;
L0=ones(Smesh,1)*0;
Nu0=1;

vec=[IL0;T0;L0;a0;R0;Nu0];

%% ODE solver
beta=beta0/BeiShu.^2;%2*kA*height*Nu0/Astar; %% could change
kn=kn0/BeiShu.^2;
[T,R]=ode45(@gene,linspace(0,tmax,Tmesh),vec);
Nu=R(:,end);
Rad=R(:,end-1);
Locations=R(:,1:Smesh);
T7=R(:,Smesh+1:2*Smesh);
Lys=R(:,2*Smesh+1:3*Smesh);
P=gamma*T7.*Lys;
AHL=R(:,3*Smesh+1);
TimeMatrix=T*ones(1,Smesh);

%%% here calculate the index equals to 0.4
A=AHL-0.4
n=2
[~,idx]=sort(abs(A(:)))
B=A(idx)
[ii,jj]=ind2sub(size(A),idx(1:n))
sortedII = sort(ii)
tpoint = sortedII(2)
%tpoint = sort(ii)(2) % the first value is the first cross point, therefore we choose the second one

tfinal = min(find(Nu-0.0001<0,1),Tmesh)+floor(27*2/BeiShu)

figure('units','normalized','outerposition',[0 0 .2 1])
subplot(1,2,1)
plot(T,AHL,'Color',[255 128 0]/255,'LineWidth',3)
xlim([0 tfinal])
ylim([0 1])
set(gca,'YTick',[0 0.4 0.8],'FontSize',15)
set(gca,'XTick',[0 50 100 150],'FontSize',15)
xlabel('Time','FontSize',15)
ylabel('AHL','FontSize',15)

subplot(2,2,2)
plot(Locations(tpoint,:)*sqrt(kn0),Lys(tpoint,:).*[ones(Smesh-1,1);0]','r','LineWidth',3)
xlim([0 2*BeiShu])
ylim([0 .4])
set(gca,'YTick',[0 .2 .4],'FontSize',15)
set(gca,'XTick',[0 1 2 3 4],'FontSize',15)
xlabel('Position','FontSize',15)
ylabel('Lysozyme','FontSize',15)

subplot(2,2,4)
plot(Locations(tfinal,:)*sqrt(kn0),Lys(tfinal,:).*[ones(Smesh-1,1);0]','r','LineWidth',3)
xlim([0 2*BeiShu])
ylim([0 .4])
set(gca,'YTick',[0 .2 .4],'FontSize',15)
set(gca,'XTick',[0 1 2 3 4],'FontSize',15)
xlabel('Position','FontSize',15)
ylabel('Lysozyme','FontSize',15)