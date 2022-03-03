function normCoeff=normPsi(psi,xmin,xmax)
normCoeff=1/sqrt(integral(matlabFunction(psi.*psi),xmin,xmax));
end
