
% range of integration
xmin = -3.5;
xmax = 3.5;
% double(int((xi(x,3)).^2,x,[xmin,xmax]))
% small number to avoid zero in denominator
eps=.1;
% Strength of repulsion (charge)
q=1;
%integration steps
steps = 25;
%initial energy
E=10;
%difference in percent with previous energy for convergence determination 
deltaEpct=1;
% Number of basis functions
K = 5;
cpar=ones(K,1)/sqrt(K);



cparArr=zeros(1,K);
EArr=zeros(1,1);
itr=0;
while deltaEpct > 0.1
    itr=itr+1
    [newE,cpar]=SCF(cpar,K,xmin,xmax,steps,eps,q);
    deltaEpct=abs(E-newE)/E.*100
    E=newE;
    cparArr(itr,:)=cpar;
    EArr(itr,1)=E;
end

Plot(K,cpar,cparArr,EArr,xmin,xmax,steps)




