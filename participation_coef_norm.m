function [PC_norm,PC_residual,between_mod_k] = participation_coef_norm(W,n_iter,par_comp,s,self,debug,verbose)

tic
%   Inputs:     W,                      binary and undirected connectivity matrix
%               Ci,                     community affiliation vector
%               n_iter (optional),      number of matrix randomizations (default = 100)
%               par_comp (optional),    0 = don't compute matrix randomizations in a parallel loop (default = 0)
%                                       1 = compute matrix randomizations in a parallel loop
%
%   Output:     PC_norm,                Normalized Participation Coefficient using randomizations
%               PC_residual,            Residual Participation Coefficient from linear regression
%               PC,                     Original Participation Coefficient
%               between_mod_k,          Between-module degree


if nargin < 3 
    n_iter = 100;
end
 
if nargin < 4 
    par_comp = 0;
end
 
n = length(W); 
Ko = sum(W,2); 
Ci = cluster_jl(W,s,self,debug,verbose);
Gc = (W~=0)*diag(Ci); 
 
Kc2 = zeros(n,1);
within_mod_k = zeros(n,1);
for i=1:max(Ci)
    Kc2 = Kc2+(sum(W.*(Gc==i),2).^2); 
    within_mod_k(Ci==i) = sum(W(Ci==i,Ci==i),2); 
end

between_mod_k = Ko - within_mod_k; 
PC = ones(n,1)-Kc2./(Ko.^2); 
PC(~Ko) = 0; % PC = 0 for nodes with no (out)neighbors

Kc2_rnd = zeros(n,n_iter); 
 
if par_comp == 0 
    reverseStr = '';
    for ii = 1:n_iter 
        ii
        if issymmetric(W) == 1 
            W_rnd = null_model_und_sign(W,5); 
        else
            W_rnd = randmio_dir(W,5); 
        end
         
        Gc_rnd = (W_rnd~=0)*diag(Ci); 
        Kc2_rnd_loop = zeros(n,1); 
 
        for iii = 1:max(Ci) 
            Kc2_rnd_loop = Kc2_rnd_loop+(((sum(W.*(Gc==iii),2)./Ko)-(sum(W_rnd.*(Gc_rnd==iii),2)./Ko)).^2);
            
        end
         
        Kc2_rnd(:,ii) = sqrt(0.5.*Kc2_rnd_loop); 
    end
else 
    parfor ii = 1:n_iter 
       if issymmetric(W) == 1 
            W_rnd = null_model_und_sign(W,5); 
        else
            W_rnd = randmio_dir(W,5); 
        end
         
        Gc_rnd = (W_rnd~=0)*diag(Ci); 
 
        for iii = 1:max(Ci) 
            Kc2_rnd_loop = Kc2_rnd_loop+(((sum(W.*(Gc==iii),2)./Ko)-(sum(W_rnd.*(Gc_rnd==iii),2)./Ko)).^2);
 
        end
         
        Kc2_rnd(:,ii) = sqrt(0.5.*Kc2_rnd_loop); 
    end
end
 
PC_norm = ones(n,1)-median(Kc2_rnd,2); 
PC_norm(~Ko) = 0; 

p = polyfit(sum(Ci==Ci')',PC,1); 
yfit = polyval(p,sum(Ci==Ci')'); 
PC_residual = PC-yfit; 
PC_residual(~Ko) = 0; 

toc

end