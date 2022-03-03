function [E,cpar]=SCF(cpar,K,xmin,xmax,steps,eps,q)
format short

syms x v xi(x,v) c(v) V(x1,x2) h(psi1) Psi

c(v)=(2/(2.*v-1));
xi(x,v)=(2.*c(v)/pi).^(1/4).*exp(-c(v).*x.^2);
h(x,v)=(-1/2).*(diff(diff(xi(x,v))))+1/2.*(x.^2).*xi(x,v);

Fh = zeros(K,K);
Fj = zeros(K,K);
F = zeros(K,K);
S = zeros(K,K);

Psi=0;
for i=1:K
    Psi=Psi+cpar(i).*xi(x,i);
end

xarr=[xmin:(xmax-xmin)/steps:xmax];
if q~=0
Jarr=intJ(Psi,xmin,xmax,eps,q,steps);
else
Jarr=zeros(1,steps+1);
end


for i=1:K
    for j=1:K
        
        hfcn=matlabFunction(xi(x,i).*h(x,j));
        Hh=trapz(xarr,hfcn(xarr));
%         if(i==1)
%             if(j==1)
%                 hfcn
%                 %xarr
%                 %fplot(hfcn, [xmin,xmax])
%                 temph=h(x,i)
%                 temphalt=halt(xi(x,i),i)
%                 fplot(matlabFunction(xi(x,i).*h(x,j)),[xmin,xmax])
%                 h(xi(x,i))
%             end
%         end


        xufcn=matlabFunction(xi(x,i));
        xvfcn=matlabFunction(xi(x,j));
        Jxiarr=xvfcn(xarr).*xufcn(xarr).*Jarr;
        HJ=trapz(xarr,Jxiarr);

        Fh(i,j)=Hh;
        Fj(i,j)=HJ;

        sfcn=matlabFunction(xi(x,i).*xi(x,j));
        S(i,j)=trapz(xarr,sfcn(xarr));

    end
end
F=Fh+Fj;

[U,A]=eig(inv(S)*F);
[eigeneps,eigenepsindex]=sort(diag(A));
minindex=eigenepsindex(1);
minvec=U(:,minindex)
Psi=0;
for i=1:K
    Psi=Psi+minvec(i).*xi(x,i);
end
normCoeff=normPsi(Psi,xmin,xmax)
cpar=minvec.*normCoeff
E=energy(K,cpar,xmin,xmax,steps,eps,q)
end
