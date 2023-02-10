function [MC,K_S,M_V,R1,R2] = MCentrality(A)

tic
K_S=kShellDecom(A); 
K_S=K_S';
len=length(A);
for i=1:len
    num_nei(i)=length(find(A(i,:)~=0));
    for j=1:len
        if A(i,j)>0
            D(i,j)=sum(A(:,j));
        else
            D(i,j)=0;
        end
    end
    
    M_V(i)=abs(sum(A(i,:)).*num_nei(i)^2/sum(D(i,:))-num_nei(i));   
    M_V(find(isnan(M_V)==1))=0;
end


for i=1:len
    R1(i)=K_S(i)/sum(K_S);
    R2(i)=M_V(i)/sum(M_V);
end

for i=1:len
    if R1(i)==0
        E_R1(i)=0;
    else
    E_R1(i)=R1(i)*log(R1(i));
    E1=-(log(len))^(-1)*sum(E_R1);
    end
    if R2(i)==0
        E_R2(i)=0;
    else
    E_R2(i)=R2(i)*log(R2(i));
    E2=-(log(len))^(-1)*sum(E_R2);
    end
end

a=(1-E1)/(2-(E1+E2));
b=1-a;

for i=1:len
    MC(i)=a*K_S(i)+b*M_V(i);
end
MC=MC';
toc
end

