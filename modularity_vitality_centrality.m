function [MVC] = modularity_vitality_centrality(A,s,self,debug,verbose)
tic
Q=modularity(A,s,self,debug,verbose);
len=length(A);
for i=1:len
    i
    Ar=A;
    Ar(i,:)=[];
    Ar(:,i)=[];
    Qr(i)=modularity(Ar,s,self,debug,verbose);
    MVC(i)=Q-Qr(i);
    
end
MVC=MVC';
toc
end
    







