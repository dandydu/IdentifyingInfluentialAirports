function [ S ] = getMatrix( str, directed )
% if it is directed network, 'directed' = 1
% if it is undirected network, 'directed' = 0
% GETMATRIX Summary of this function goes here
% Detailed explanation goes here
if nargin == 1
    directed = 1;
end

D = textread(str);
l = length(D);
n = max(max(D(:,1:2)));
t = min(min(D(:,1:2)));
[~,m] = size(D);
if m == 2
    D = [D,ones(l,1)];
end
if t == 0
    t = 1;
else
    t = 0;
end
S = zeros(n+t);
for i = 1:length(D)
    S(D(i,1)+t,D(i,2)+t) = D(i,3);
end
if directed == 0
    S = S+S';
end


end

