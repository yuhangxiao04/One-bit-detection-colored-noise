function p=orthant(C,np)

%The function to compute the orthant probability of a correlation matrix C.

%C is the correlation matrix of the noise.
%np is the number of repeated times in computing the orthant probabilities.



[m,~]=size(C);

u=zeros(1,m);

thre=zeros(1,m);

p_mont=zeros(1,np);

for i=1:np
 
    p_mont(i)=mvncdf(thre,u,C);

end

p=mean(p_mont);