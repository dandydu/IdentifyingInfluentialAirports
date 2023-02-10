function [deg] = degrees_und(CIJ)

CIJ = double(CIJ~=0);

deg = sum(CIJ);

