function plot(K,cpar,cparArr,EArr,xmin,xmax,steps)
%K=4;
%cpar=[0.0303,-0.6785,1.7564,-0.1202];%result for q=1
%cpar=[-0.0105,0.3462,0.8554,-0.1914];%result for q=0 
%xmin=-3.5;
%xmax=3.5;
%steps=25;

syms xi(x,v) c(v) V(x1,x2) h(psi1) Psi1p

c(v)=(2/(2.*v-1));
xi(x,v)=(2.*c(v)/pi).^(1/4).*exp(-c(v).*x.^2);

Psi=0;
for i=1:K
    Psi=Psi+cpar(i).*xi(x,i);
end
%check normalization
normCoeffSCF=normPsi(Psi,xmin,xmax)
Psi2fcn=matlabFunction(Psi.*Psi);
xarr=[xmin:(xmax-xmin)/steps:xmax];
figure();
fplot(Psi2fcn,[xmin,xmax])
xlabel("x(a.u)");
ylabel("|\Psi(x)|^2");

hold on
Psi1p=(1/pi).^(1/4).*exp(-x.^2/2);
normCoeffGRD=normPsi(Psi1p,xmin,xmax)
Psi1p2fcn=matlabFunction(Psi1p.^2);
fplot(Psi1p2fcn,[xmin,xmax])
title('','Interpreter','latex');
label1 = '$ |\Psi_{SCF}|^2$';
label2 = '$ |\Psi_{Ground}|^2$';
legend(label1,label2,'Interpreter','latex')
hold off

figure();
Psifcn=matlabFunction(Psi);
Psi1pfcn=matlabFunction(Psi1p);
fplot(Psifcn,[xmin,xmax])
xlabel("x(a.u)");
ylabel("|\Psi(x)|^2");
hold on
Psi1pfcn=matlabFunction(Psi1p);
fplot(Psi1pfcn,[xmin,xmax])
title('','Interpreter','latex');
label1 = '$ \Psi_{SCF}$';
label2 = '$ \Psi_{Ground}$';
legend(label1,label2,'Interpreter','latex')
hold off

figure();
itrs=size(cparArr,1);
xitr=[1:1:itrs];
labels=strings(K);
for i = 1:K
    labels(i)=sprintf('$C_{%d}$', i);
    plot(xitr,cparArr(:,i));
    if i==1
        xlabel("Iteration");
        hold on
    end
end
legend(labels,'Interpreter','latex');

figure();
plot(xitr,EArr(:,1));
xlabel("Iteration");
ylabel("Energy (a.u.)")

end








