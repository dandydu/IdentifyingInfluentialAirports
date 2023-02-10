function [F]=F_SI(net_adj,net_edge)

col=length(net_adj);
alpha = 1;
tt = 1; % Record the number of steps
maxStep = 10;
Sta = 1000; 
Trace = zeros(maxStep,Sta);
Node_Inf = zeros(Sta,1);

Fi=zeros(10,Sta);
for i=1:col
    IOi=[i]
   
    for k = 1:Sta
        [ F] = SI(net_edge, alpha, IOi, tt*maxStep);
        Fi(i,k)=F(maxStep);
    end
    
end

F=sum(Fi,2)./Sta;
end