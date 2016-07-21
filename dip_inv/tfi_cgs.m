function chi = tfi_cgs(tfs, mask, vox, z_prjs, L2_sus_weight, L2_lfs_weight)

% total field inversion
% ||M * F_{-1} * D * F * chi_t - M * B_t||_2 
%   + L2_sus_weight*||M * chi_t||_2
%	  + L2_lfs_weight * ||F_{-1} * D * D * M *chi_t||_2
%   

[Nx,Ny,Nz] = size(tfs);
imsize = size(tfs);

% create K-space filter kernel D
%%%%% make this a seperate function in the future
FOV = vox.*[Nx,Ny,Nz];
FOVx = FOV(1);
FOVy = FOV(2);
FOVz = FOV(3);

x = -Nx/2:Nx/2-1;
y = -Ny/2:Ny/2-1;
z = -Nz/2:Nz/2-1;
[kx,ky,kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
% D = 1/3 - kz.^2./(kx.^2 + ky.^2 + kz.^2);
D = 1/3 - (kx.*z_prjs(1)+ky.*z_prjs(2)+kz.*z_prjs(3)).^2./(kx.^2 + ky.^2 + kz.^2);
D(floor(Nx/2+1),floor(Ny/2+1),floor(Nz/2+1)) = 0;
D = fftshift(D);



%%%%%%%%%%%%%%% cgs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = ifftn(conj(D).*fftn(mask.*tfs));
b = b(:);

chi = cgs(@Afun, b, 1e-6, 500);

chi = reshape(chi,imsize);



% nested function
function y = Afun(x)
    % y = H'*(H*x) + tik_reg*x;
    x = reshape(x,imsize);
    y = ifftn(D.*fftn(mask.*ifftn(D.*fftn(x)))) + L2_sus_weight*mask.*x + L2_lfs_weight*mask.*(ifftn(D.*D.*fftn(mask.*x)));

    y = y(:);
end

end