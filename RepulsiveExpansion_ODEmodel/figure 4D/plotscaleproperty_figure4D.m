clear all
close all
global Nstar Kphi Kphi2 Smesh gamma beta Tstar Pstar da dT dL sigma0 kT kL m expphi expphi2 kn

%% Parameter Statement

L=1;  % length of interval
tmax=251; % time integration
Tmesh=200;
Smesh=100;
N=400+1;   % number of grid points; spatial step size h=L/(N-1)

%%%%%%%% the parameters %%%%%%%%%%%%%%%%%
kT = 8500;
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
a0=0;   %change
Nstar=0.9;
Kphi=0.1;
expphi=2;
Kphi2=0.1;
expphi2=4;
m = 4;
beta0=81;
kn0=25;
tmesh=160;

BeiShuD=[.4 .6 .8 1 1.2 1.4 1.6 1.8 2];
ColonyRadius=zeros(1,size(BeiShuD,2));
RingWidth=zeros(1,size(BeiShuD,2));

for i=1:size(BeiShuD,2)
    BeiShu=BeiShuD(i);
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
kn=kn0/BeiShu.^2;

[T,R]=ode45(@gene_20171012,linspace(0,tmax,Tmesh),vec);
Nu=R(:,end);
Rad=R(:,end-1);
Locations=R(:,1:Smesh);
T7=R(:,Smesh+1:2*Smesh);
Lys=R(:,2*Smesh+1:3*Smesh);
P=gamma*T7.*Lys;
AHL=R(:,3*Smesh+1);
TimeMatrix=T*ones(1,Smesh);

tfinal=min(find(Nu-0.0001<0,1),Tmesh)+floor(27*2/BeiShu);
if size(tfinal,1)==0
  tfinal=tmesh
else
  end

mid=find(Lys(tfinal,:)==min(Lys(tfinal,:)));
if mid==1
  RW=0;
else
  RW=(Rad(end)-Locations(tfinal,mid))*sqrt(kn0);
  end
  
ColonyRadius(i)=Rad(end)*sqrt(kn0);
RingWidth(i)=RW;
end

range=2:8;
p1 = polyfit(BeiShuD(range),ColonyRadius(range),1);
slope1=p1(1,1);
intercept1=p1(1,2);
p2 = polyfit(BeiShuD(range),RingWidth(range),1);
slope2=p2(1,1);
intercept2=p2(1,2);
z=0:0.1:4;
plot(z,slope1.*z+intercept1,'Color',[0 128 0]/255,'LineWidth',3)
hold on
plot(BeiShuD(1:end),ColonyRadius(1:end),'o','Color',[0 128 0]/255,'MarkerFaceColor','w','MarkerSize',15,'LineWidth',3)
hold on
plot(BeiShuD(range),ColonyRadius(range),'o','Color',[0 128 0]/255,'MarkerFaceColor',[0 128 0]/255,'MarkerSize',15,'LineWidth',3)
hold on
plot(z,slope2.*z+intercept2,'r','LineWidth',3)
hold on
plot(BeiShuD(1:end),RingWidth(1:end),'or','MarkerFaceColor','w','MarkerSize',15,'LineWidth',3)
hold on
plot(BeiShuD(range),RingWidth(range),'or','MarkerFaceColor','r','MarkerSize',15,'LineWidth',3)

xlim([0 4]);
ylim([0 5]);
set(gca,'YTick',0:2:5,'FontSize',20)
set(gca,'XTick',0:1:4,'FontSize',20)
xlabel('Domain radius','FontSize',20)
ylabel('Distance','FontSize',20)
box on
box square