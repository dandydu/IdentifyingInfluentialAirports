function [CC]=closeness_centrality(Net,Net2)

%INPUT��Net is adjaceny matrix, Net2 is adjaceny list
%OUTPUT��the closeness of node

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


