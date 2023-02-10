function [id,od,deg] = degrees_dir(CIJ)

CIJ = double(CIJ~=0);

% compute degrees
id = sum(CIJ,1); 
od = sum(CIJ,2)'; 
deg = id+od; 