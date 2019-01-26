function L=cholesky(A)  
n=size(A); 
L=zeros(n); 
L(1,1)=sqrt(A(1,1)); 
for i=2:1:n 
    L(i,1)=A(i,1)/L(1,1); 
end 
for j=2:1:n 
    sum1=0; 
    for k=1:1:j-1 
        sum1=sum1+L(j,k)^2; 
    end 
    L(j,j)=sqrt(A(j,j)-sum1); 
    for i=j+1:1:n; 
        sum2=0; 
        for k=1:1:j-1 
            sum2=sum2+L(i,k)*L(j,k); 
        end 
        L(i,j)=(A(i,j)-sum2)/L(j,j); 
    end 
end 