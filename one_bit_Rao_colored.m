function [T,s1,s2] = one_bit_Rao_colored(A,B,Y,L,p,Var)

%The function to construct the detector according to the sample 
%and the parameters computed according to the noise covariance matrix; 


    [m,n]=size(Y);

    Y1=[real(Y);imag(Y)];

    u=2.^fliplr(0:2*m-1);

    z=u*(Y1+1)/2+1;

    LP_inv=L*diag(p)^(-1);

    d=zeros(2*m,n);

    for j=1:n      
        d(:,j) = LP_inv(:,z(j));
    
    end

    s1=sum(sum(d.*A));

    s2=sum(sum(d.*B));


    T = (s1^2+s2^2)/Var;


