function tau = KendallTau(lista, listb)
%Kendall Tau
lista=double(full(lista));
listb=double(full(listb));
n=length(lista);

num=0;
for ii=1:n-1,
    for jj=ii+1:n,
        num=num+sign(lista(ii)-lista(jj))*sign(listb(ii)-listb(jj));
        %fprintf('%d %d\n',ii,jj);
    end
end
tau=2*num/(n^2-n);

end

