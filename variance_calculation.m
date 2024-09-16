function Var = variance_calculation(A,B,L,p)

%The function to calculate theoretical variance.

Var1=trace((A'*L)'*(A'*L)*diag(p)^(-1));

Var2=trace((B'*L)'*(B'*L)*diag(p)^(-1));

Var=(Var1+Var2)/2;