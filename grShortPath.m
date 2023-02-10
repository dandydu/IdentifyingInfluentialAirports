function [dSP,sp]=grShortPath(E,s,t)

if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

dSP=ones(n)*inf; % initial distances
dSP((E(:,2)-1)*n+E(:,1))=E(:,3);

for j=1:n,
  i=setdiff((1:n),j);
  dSP(i,i)=min(dSP(i,i),repmat(dSP(i,j),1,n-1)+repmat(dSP(j,i),n-1,1));
end
sp=[];
if (nargin<3)|(isempty(s))|(isempty(t)),
  return
end
s=s(1);
t=t(1);
if (~(s==round(s)))|(~(t==round(t)))|(s<1)|(s>n)|(t<1)|(t>n),
  error(['s and t must be integer from 1 to ' num2str(n)])
end
if isinf(dSP(s,t)), % t is not accessible from s
  return
end
dSP1=dSP;
dSP1(1:n+1:n^2)=0;
l=ones(m,1);
sp=t;
while ~(sp(1)==s),
  nv=find((E(:,2)==sp(1))&l);
  vnv=abs((dSP1(s,sp(1))-dSP1(s,E(nv,1)))'-E(nv,3))<eps*1e3;
  l(nv(~vnv))=0;
  if all(~vnv),
    l(find((E(:,1)==sp(1))&(E(:,2)==sp(2))))=0; 
    sp=sp(2:end);
  else
    nv=nv(vnv);
    sp=[E(nv(1),1) sp];
  end
end
return