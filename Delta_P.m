function P_delta = Delta_P(H_tilde,C_delta)

%The function to calculate the change of P with regard to the change in C,
%with the aid of H_tild.

v=zeros(1,2*m*(m-1));
kk=1;
for j=1:2*m-1
    for s=j+1:2*m
        if s~=j+m
            v(kk)=C_delta(j,s);
            kk=kk+1;
        end
    end
end


P_delta=diag(v*H_tilde);