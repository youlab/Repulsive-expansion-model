function out=gene_20171012(t,vec)
global Nstar Kphi Kphi2 Smesh gamma beta Tstar Pstar da dT dL sigma0 kT kL m expphi expphi2 kn
%% unpack the vectors used in the function

radius=vec(1:Smesh);
T=vec(Smesh+1:2*Smesh);
Lys=vec(2*Smesh+1:3*Smesh);
P=gamma*T.*Lys;
a=vec(3*Smesh+1);
R=vec(3*Smesh+2);
Nu=vec(3*Smesh+3);

NTerm=(Nu)./(Nu+Nstar);
GK=Kphi.^expphi./(Kphi.^expphi+(max(R-radius,0)).^expphi);
GT=Kphi2.^expphi2./(Kphi2.^expphi2+(radius).^expphi2);%growth term
phi=NTerm.*GT;%./(1+0.1*T);% growth/divide rate
r_mat=radius.*phi;
F_R=1/R*trapz(radius,r_mat);
F_r(1,1)=0;
for i=2:Smesh
    F_r(i,1)=1./radius(i)*trapz(radius(1:i),r_mat(1:i));
end
alpha=T./(1+T)./(1+Tstar/Pstar*P).*GK;
omega=T./(1+T).*(a.^m./(1+a.^m)).*GK;

a_mat=alpha.*radius;

%% T %%%%%%%%%%%%
F_T=1./(1+gamma*Lys+gamma*T).*((1+gamma*T).*(kT/Tstar/sigma0*alpha-dT/sigma0*T)-gamma*T.*(kL/Tstar/sigma0*omega-dL/sigma0*Lys));

%% Lys %%%%%%%%%%%
F_L=1./(1+gamma*Lys+gamma*T).*((1+gamma*Lys).*(kL/Tstar/sigma0*omega-dL/sigma0*Lys)-gamma*Lys.*(kT/Tstar/sigma0*alpha-dT/sigma0*T));

%% AHL concentration %%%%%%%%%%%
F_a=beta.*trapz(radius,a_mat)-da/sigma0.*a; %% one cell is 1um=1e-4cm

%% nutrient concentration %%%%%%%%%%%
F_Nu=-kn.*trapz(radius,r_mat);

%% Assemble them into the output vector
out=[F_r;F_T;F_L;F_a;F_R;F_Nu];

