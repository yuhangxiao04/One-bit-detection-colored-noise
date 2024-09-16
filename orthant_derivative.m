function L=orthant_derivative(C,np)

%The function to compute all possible values of the columns of the matrix "L" in the paper. 

%C is the correlation matrix of the noise.
%np is the number of repeated times in computing the orthant probabilities.

[m,~]=size(C);

m=m/2;

L=zeros(2*m,2^(2*m));

T=[zeros(m,m),-eye(m);
   eye(m),zeros(m,m)];

T1=[zeros(m,m),eye(m);
   eye(m),zeros(m,m)];

for k=1:2^(2*m)
   
    cc=str2num(strread(dec2bin(k-1,2*m), '%c'));

    cc=cc*2-1;

    cc1=T*cc;

    cc2=T*cc1;

    cc3=T*cc2;

    k1=2.^fliplr(0:2*m-1)*(cc1+1)/2+1;

    k2=2.^fliplr(0:2*m-1)*(cc2+1)/2+1;

    k3=2.^fliplr(0:2*m-1)*(cc3+1)/2+1;

    if k==min([k,k1,k2,k3])

    parfor i=1:2*m

        L=diag(cc)*C*diag(cc);

        t=L(i,:);

        t(i)=[];

        L(i,:)=[];

        L(:,i)=[];

        F=L-(t'*t)/C(i,i);

        L(i,k) = orthant(F,np);

    end

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

    if L(1,k)~=0

        L(:,k1) =  T1*L(:,k);

        L(:,k2) =  L(:,k);

        L(:,k3) =  T1*L(:,k);

    end
end

L=L/sqrt(2*pi);

for j=1:2^(2*m)

    cc=str2num(strread(dec2bin(j-1,2*m), '%c'));

    cc=cc*2-1;

    L(:,j)=diag(cc)*L(:,j);

end