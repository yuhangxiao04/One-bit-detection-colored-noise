function H_tilde = H_tilde_fast(C,np)

%A support function to calculate the orthant probabilities in the mismtached case. 



[m,~]=size(C);

m=m/2;

T=[zeros(m,m),-eye(m);
   eye(m),zeros(m,m)];
dd=(m-1)*m/2;

H_tilde=zeros(2*m*(m-1),2^(2*m-1));
L=zeros(dd,2^(2*m-1));


kk=1;
kk1=1;
for j=1:m-1

    for k=j+1:m        

            parfor i=1:2^(2*m-1)
                cc=str2num(strread(dec2bin(i-1,2*m), '%c'));

                cc=cc*2-1;

                C_temp=diag(cc)*C*diag(cc);

                L(kk1,i) = orthant_de(C_temp,j,k,np);
                H_tilde(kk,i) = cc(j)*cc(k)*L(kk1,i);
            end

           kk=kk+1;
           kk1=kk1+1;
    end
    kk=kk+(m-1);

end

H_tilde=[H_tilde,fliplr(H_tilde)];
L=[L,fliplr(L)];

kk=dd*3+1;

for j=1:m-1

    for k=j+1:m

        for i=1:2^(2*m-1)

            cc=str2num(strread(dec2bin(i-1,2*m), '%c'));

            cc=cc*2-1;

            cc1=T*cc;

            k1=2.^fliplr(0:2*m-1)*(cc1+1)/2+1;

            H_tilde(kk,i) = cc(j+m)*cc(k+m)*L(kk-dd*3,k1);

            H_tilde(kk,2^(2*m)-i+1) = H_tilde(kk,i);

        end

        kk=kk+1;

    end

end



L1=zeros(dd,2^(2*m-1));
kk1=1;

for j=1:m-1

    for k=j+1+m:2*m        

            parfor i=1:2^(2*m-1)
                cc=str2num(strread(dec2bin(i-1,2*m), '%c'));

                cc=cc*2-1;

                C_temp=diag(cc)*C*diag(cc);

                L1(kk1,i) = orthant_de(C_temp,j,k,np);
            end

           kk1=kk1+1;

    end

end

L1=[L1,fliplr(L1)];

O1=zeros(m-1,m-1);

ll=1;
for i=1:m-1
    for j=i:m-1
        O1(i,j)=ll;
        ll=ll+1;
    end
end

O2=O1';

ind=zeros(1,m*(m-1));

ll=1;
for i=1:m-1
    for j=1:i
        ind(ll)=O2(i,j);
        ll=ll+1;
    end
end


kk=m;
kk1=1;
kk2=1;

for j=1:m

    for k=m+1:2*m

        if k<j+m

            for i=1:2^(2*m-1)

                cc=str2num(strread(dec2bin(i-1,2*m), '%c'));

                cc=cc*2-1;

                cc1=T*cc;

                k1=2.^fliplr(0:2*m-1)*(cc1+1)/2+1;

                H_tilde(kk,i) = cc(j)*cc(k)*L1(ind(kk1),k1);

                H_tilde(kk,2^(2*m)-i+1) = H_tilde(kk,i);

            end

            kk1=kk1+1;

            kk=kk+1;

        end

        if k>j+m


            for i=1:2^(2*m-1)

                cc=str2num(strread(dec2bin(i-1,2*m), '%c'));

                cc=cc*2-1;

                H_tilde(kk,i) = cc(j)*cc(k)*L1(kk2,i);

                H_tilde(kk,2^(2*m)-i+1) = H_tilde(kk,i);

            end

            kk2=kk2+1;

            kk=kk+1;

        end

    end

        kk=kk+m-j-1;

end




