function d = simple_dijkstra(adj,s)

n=length(adj);
d = inf*ones(1,n);
d(s) = 0;
T = 1:n; 

while not(isempty(T))
    [dmin,ind] = min(d(T));
    for j=1:length(T)
        if adj(T(ind),T(j))>0 & d(T(j))>d(T(ind))+adj(T(ind),T(j))
            d(T(j))=d(T(ind))+adj(T(ind),T(j));
        end
    end 
    T = setdiff(T,T(ind));
end