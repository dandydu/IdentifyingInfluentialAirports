function [DCC] = degree_clustering_coefficient(A)

tic
len=length(A);
for i=1:len
    for j=1:len
        if A(i,j)>0
            D(i,j)=sum(A(:,j));
        else
            D(i,j)=0;
        end
    end
    I_D(i)=sum(A(i,:))+sum(D(i,:));
end


for i=1:len
    if issymmetric(A)==1 
        N(:,i)=[neighbors(graph(A),i);zeros(len-length(neighbors(graph(A),i)),1)];
        
    else
        N(:,i)=[predecessors(digraph(A),i);zeros(len-length(predecessors(digraph(A),i)),1)];
        
    end
end

for i=1:len
    eval(['N' num2str(i) '=N(:,N(find(N(:,i)~=0),i));']);
end

NN={};

for j=1:len
    NN{j}=eval(['N' num2str(j)]);
end

for i=1:len
NN{1,i}(find(NN{1,i}==i))=[];
NN{1,i}(find(NN{1,i}==0))=[];
NN{1,i}=unique(NN{1,i});
end

ccfs = clustering_coefficients(sparse(A)); 

for i=1:len
    for j=1:length(NN{1,i})
        CCFSS(i,j)=ccfs(NN{1,i}(j));
        I_C(i)=exp(-ccfs(i))*sum(CCFSS(i,:));
        
    end
end

L1=length(I_C);
L2=length(NN);
LL=L1-L2;
if LL<0
    I_C(L1+1:L2)=0;
end

for i=1:len
    R_D(i)=I_D(i)/sqrt(sum(I_D.*I_D,2));
    R_C(i)=I_C(i)/sqrt(sum(I_C.*I_C,2));
end

for i=1:len
    if R_D(i)==0
        E_D(i)=0;
    else
    E_D(i)=R_D(i)*log(R_D(i));
    E1=-(log(len))^(-1)*sum(E_D);
    end
    if R_C(i)==0
        E_C(i)=0;
    else
    E_C(i)=R_C(i)*log(R_C(i));
    E2=-(log(len))^(-1)*sum(E_C);
    end
end

a=(1-E1)/(2-(E1+E2));
b=1-a;

for i=1:len
    
    DCC(i)=a*I_D(i)+b*I_C(i);
    DCC=DCC';
end
toc
end

