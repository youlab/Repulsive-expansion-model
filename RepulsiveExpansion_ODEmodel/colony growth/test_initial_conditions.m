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

% %% test 1: initial colony size
% figure('units','normalized','outerposition',[0 0 .5 .8])
% for N=[10 50 100 200 400]
%     %% Define Initial Conditions
%     R0=L/(N-1);
%     IL=1:1:Smesh;
%     IL0=(exp(log(R0+1)/(Smesh-1)*(IL-1))-1)';
% 
%     Nu0=1;
%     vec=[IL0;R0;Nu0];
%  
%     %% ODE solver
%     [T,R]=ode45(@gene,linspace(0,tmax,Tmesh),vec);
%     Nu=R(:,end);
%     Rad=R(:,end-1);
%     Locations=R(:,1:Smesh);
%     TimeMatrix=T*ones(1,Smesh);
%     subplot(4,3,1:3)
%     plot([IL0; IL0(end)],[ones(Smesh,1)*1; 0],'LineWidth',3)
%     hold on
%     xlim([0 0.12])
%     ylim([0 1.05])
%     set(gca,'XTick',[0 0.1 .2 .3 .4 .5],'FontSize',20)
%     set(gca,'YTick',[0 1],'FontSize',20)
% %     xlabel('Position','FontSize',15)
% %     ylabel('Cell Density','FontSize',15)
%     subplot(4,3,4:12)
%     plot(T,Rad,'LineWidth',3)
%     hold on
%     xlim([0 tmax+1])
%     ylim([0 1.6])
%     set(gca,'XTick',[0 10 20],'FontSize',20)
%     set(gca,'YTick',[0 0.5 1 1.5 2],'FontSize',20)
%     xlabel('Time','FontSize',20)
% %     ylabel('Colony Radius','FontSize',15)  
% %     set(gca,'DataAspectRatio',[1 7 1])
%     box off
% end



%% test 2: initial position selection
figure('units','normalized','outerposition',[0 0 1 .25])
ha= tight_subplot(1,4,[.01 .01],[.45 .25],[.05 .01]); %(Nh, Nw, gap, marg_h, marg_w)
% ha2= tight_subplot(1,4,[.01 .01],[.25 .1],[.05 .01]);
for i = [1 2 3 4]
    % Define Initial Conditions
    R0=L/(N-1);
    if i==1
        IL0=linspace(0,R0,Smesh)';
    elseif i==2
        IL=linspace(0,R0,Smesh)';
        IL0=2*IL-IL.^2/R0;
    elseif i==3
        IL=linspace(0,R0,Smesh)';
        IL0=IL.^2/R0;
    else
        arbi=3;
        I=linspace(log(R0)/log(arbi),-20,Smesh-2)';
        IL0=[0; R0-arbi.^(I); R0];
    end
        
    Nu0=1;
    vec=[IL0;R0;Nu0];
    
%     % ODE solver
%     [T,R]=ode45(@gene,linspace(0,tmax,Tmesh),vec);
%     Nu=R(:,end);
%     Rad=R(:,end-1);
%     Locations=R(:,1:Smesh);
%     TimeMatrix=T*ones(1,Smesh);

    axes(ha(i))
    plot([IL0; IL0(end)],[ones(Smesh,1)*1; 0],'r.','MarkerSize',5)
    xlim([0 0.0101])
    ylim([0 1.5])
    if i==1
    set(gca,'XTick',[0 0.005 0.01],'FontSize',15)
    set(gca,'YTick',[0 1],'FontSize',15)
    xlabel('Position','FontSize',15)
    ylabel('Cell Density','FontSize',15)
    else
        set(gca,'YTick',[])
        set(gca,'XTick',[])
    end

%     axes(ha2(i))
%     for ii=1:5:Smesh
%         plot(R(:,ii),T,'LineWidth',1)
%         hold on
%     end
%     plot(Rad,T,'k','LineWidth',2)
%     hold on
%     xlim([0 1.5])
%     ylim([0 tmax+1])
%     if i==1
%         set(gca,'YTick',[0 10 20],'FontSize',15)
%         set(gca,'XTick',[0 0.5 1 1.5],'FontSize',15)
%         xlabel('Position','FontSize',15)
%         ylabel('Time','FontSize',15)
%     else
%         set(gca,'YTick',[])
%         set(gca,'XTick',[])
%     end
%     axis square
end
