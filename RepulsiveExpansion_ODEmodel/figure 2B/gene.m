function out=gene(t,vec)
global Nstar  Kphi2 Smesh expphi2 kn
%% unpack the vectors used in the function

radius=vec(1:Smesh);
R=vec(Smesh+1);
Nu=vec(Smesh+2);

NTerm=(Nu)./(Nu+Nstar);
GT=Kphi2.^expphi2./(Kphi2.^expphi2+(R-radius).^expphi2);%growth term
phi=NTerm.*GT; % growth/divide rate
r_mat=radius.*phi;
F_R=1/R*trapz(radius,r_mat);
F_r(1,1)=0;
for i=2:Smesh
    F_r(i,1)=1./radius(i)*trapz(radius(1:i),r_mat(1:i));
end
%% nutrient concentration %%%%%%%%%%%
F_Nu=-kn.*trapz(radius,r_mat);

%% Assemble them into the output vector
out=[F_r;F_R;F_Nu];

