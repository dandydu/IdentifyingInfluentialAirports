function [Result]=node_TOPSIS(network)

tic

[Occurence,Severity,Detection]=FMEA_RPN(network);
CN(:,1)=Occurence;
CN(:,2)=Severity;
CN(:,3)=Detection;

sum1=sqrt(sum(CN(:,1).^2));
sum2=sqrt(sum(CN(:,2).^2));
sum3=sqrt(sum(CN(:,3).^2));

col=length(CN(:,1));
for i=1:col
        B(i,1)=1/3.*CN(i,1)/sum1;
end
for i=1:col
        B(i,2)=1/3.*CN(i,2)/sum2;
end
for i=1:col
        B(i,3)=1/3.*CN(i,3)/sum3;
end

C1=B(:,1);
C2=B(:,2);
C3=B(:,3);

max1=max(B(:,1));
min1=min(B(:,1));

max2=max(B(:,2));
min2=min(B(:,2));

max3=max(B(:,3));
min3=min(B(:,3));

col=length(CN(:,1));
for i=1:col
Dmax(i)=sqrt((max1-C1(i)).^2+(max2-C2(i)).^2+(max3-C3(i)).^2);
Dmin(i)=sqrt((min1-C1(i)).^2+(min2-C2(i)).^2+(min3-C3(i)).^2);
end

Result=Dmin./(Dmin+Dmax);
Result=Result';
toc
end
