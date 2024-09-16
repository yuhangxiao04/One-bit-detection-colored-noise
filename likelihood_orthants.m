function p=likelihood_orthants(C,np)

%The function to compute all possible values of the diagonal elements of the matrix "P" in the paper. 

%C is the correlation matrix of the noise.
%np is the number of repeated times in computing the orthant probabilities.



[m,~]=size(C);

m=m/2;

T=[zeros(m,m),-eye(m);
   eye(m),zeros(m,m)];


p=zeros(1,2^(2*m));

parfor k=1:2^(2*m)

    cc=str2num(strread(dec2bin(k-1,2*m), '%c'));

    cc=cc*2-1;

    cc1=T*cc;

    cc2=T*cc1;

    cc3=T*cc2;

    k1=2.^fliplr(0:2*m-1)*(cc1+1)/2+1;

    k2=2.^fliplr(0:2*m-1)*(cc2+1)/2+1;

    k3=2.^fliplr(0:2*m-1)*(cc3+1)/2+1;

    if k==min([k,k1,k2,k3])

        p(k) = orthant(diag(cc)*C*diag(cc),np);

    end

end



for k=1:2^(2*m)

    cc=str2num(strread(dec2bin(k-1,2*m), '%c'));

    cc=cc*2-1;

    cc1=T*cc;

    cc2=T*cc1;

    cc3=T*cc2;

    k1=2.^fliplr(0:2*m-1)*(cc1+1)/2+1;

    k2=2.^fliplr(0:2*m-1)*(cc2+1)/2+1;

    k3=2.^fliplr(0:2*m-1)*(cc3+1)/2+1;

    if p(k)~=0

        p(k1) = p(k);

        p(k2) = p(k);

        p(k3) = p(k);
    end
end


