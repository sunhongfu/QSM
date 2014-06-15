% LBV.m
% Full Multigrid (FMG) solver for the Laplacian Boundary Value problem
% written by Dong Zhou (zhou.dong@gmail.com)
% this is a wrapper for mexMG.cpp
% created:   : 6.12.2013
% last modify: 1.22.2014
%
% please cite
% Background field removal by solving the Laplacian boundary value problem
% Dong Zhou, Tian Liu, Pascal Spincemaille, Yi Wang, NMR BioMed 2014
% DOI: 10.1002/nbm.3064
% 
% http://onlinelibrary.wiley.com/doi/10.1002/nbm.3064/abstract
% ----------------------------
%   [fL] = LBV(iFreq,Mask,matrix_size,voxel_size,tol,depth,peel,N1,N2,N3)
% 
%   output
%   fL - the local field
% 
%   input
%   iFreq - total magnetic field
%   Mask - ROI
%   matrix_size - dimension of the 3D image stack
%   voxel_size  - dimensions of the voxels in unit of mm
%   tol  -  iteration stopping criteria on the coarsest grid
%           try 0.01 or so
%   peel -  number of boundary layers to be peeled off
%           this is a quick-and-dirty way to get rid of 
%           bad voxels at the ROI boundary
%   depth - number of length scales. 
%           The largest length scale is 2^depth * voxel size.
%   N1 - iterations on each depth before the recursive call
%   N2 - iterations on each depth after the recursive call
%   N3 - iterations on the finest scale after the FMG is finished.
%
%

function [fL] = LBV(iFreq,Mask,matrix_size,voxel_size,tol,peel,depth,N1,N2,N3)
    if (nargin <6)
        peel = 0;
    end
    if (nargin < 7)
        depth = -1;
    end
    if (nargin < 8)
        N1 = 30;
    end
    if (nargin < 9)
        N2 = 100;
    end
    if (nargin < 10)
        N3 = 100;
    end

    % data type conversion
    matrix_size = double(matrix_size);
    voxel_size = double(voxel_size);
    n_vox = prod(matrix_size);
    fT = double(reshape(iFreq, [1,n_vox]));
    mask = double(reshape(Mask, [1,n_vox]));

    % c++ interface
    tmp = mexMGv6(fT,mask,matrix_size,voxel_size,tol,depth,peel,N1,N2,N3);
    fL = reshape(tmp,matrix_size);

end

