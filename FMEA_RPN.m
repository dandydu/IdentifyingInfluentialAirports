function [Occurence,Severity,Detection]=FMEA_RPN(network)

A=network;
L=length(A);

%   Occurence based on indegree
alph=0.5;
AA = double(A~=0);
indc = sum(AA,1);
Occurence=(indc/(L-1)).^alph;
Occurence=Occurence';


%   Severity based on effective distance closeness centrality
a=sum(A,2);
b=meshgrid(a)';
d=1-log(A./b);
d(isnan(d))=0;
d(d==inf)=0;
d=d';

col=length(d);
D=zeros(col);
C=zeros(1,col);
for i=1:col
    %i
    D(i,:) = simple_dijkstra(d,i);
end

D(D==inf)=0;

for i=1:col
    C(i)=1./sum(D(:,i));
end

C(isnan(C))=0;
C(C==inf)=0;
Severity=C';
Severity(find(Severity==1))=0;


%   Detection based on node_entropy
for i=1:L
    for j=1:L
        if A(i,j)>0
            D(i,j)=sum(A(j,:));
        else
            D(i,j)=0;
        end
    end
    E(i,:)=D(i,:)/sum(D(i,:));
end

for i=1:L
    for j=1:L
        if E(i,j)>0
            H(i,j)=sum(-E(i,j)*log2(E(i,j)));
        else
            H(i,j)=0;
        end
    end
    HH(i)=sum(H(i,:));
end
Detection=HH';
end
