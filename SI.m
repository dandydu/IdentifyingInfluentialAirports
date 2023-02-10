function [F] = SI(net,alpha,IO,maxStep)
% net, adjaceny matrix
% alpha, infection probability
% IO, infected node
% maxStep, the max step of infection

if size(net,2)==2
    net(:,3)=1;
end

n=max(max(net(:,1:2)));
m=length(net);
prob=(net(:,3)./(2*max(net(:,3)))).^alpha; 
state=zeros(n,1);
state(IO)=1;
step =1;
F=zeros(maxStep,1);
F(step)=1;
for i =1 :maxStep-1
    edgeState=state(net(:,1)); 
    state(net((edgeState.*rand(m,1)>1-prob),2))=1; 
    step=step+1;
    F(step)=sum(state);
    
end

end

