function J0int=intJ(Psi,xmin,xmax,eps,q,steps) 
    syms J xi(x,v) c(v) Psi2 V(x1,x2)
%     xmin=-100;
%     xmax=100;
%     K=4;
%     cpar=ones(4,1)/sqrt(K);
%     eps=.1;
%     q=1;

    V(x1,x2)=q.^2/(abs(x1-x2)+eps);

    normCoeff=normPsi(Psi,xmin,xmax);
%     J=int(Psi2.*V(x1,x2),x1,[xmin,xmax])
    xarr=[xmin:(xmax-xmin)/steps:xmax];
    J0=matlabFunction(normCoeff.*normCoeff.*Psi.*Psi.*V(xarr,x));
    J0int=integral(J0,xmin,xmax,'ArrayValued',true);
   
end
