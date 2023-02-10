function [CC]=closeness_centrality(Net,Net2)

%INPUT£ºNet is adjaceny matrix, Net2 is adjaceny list
%OUTPUT£ºthe closeness of node

tic
col=length(Net);
dd=grShortPath(Net2);
dd(logical(eye(size(dd))))=0;
for i=1:col
    i
    CC(i)=1./sum(dd(i,:));
end
CC=CC';
toc
end


