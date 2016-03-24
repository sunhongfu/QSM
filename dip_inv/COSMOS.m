% run QSM individually on each orientation




% define common space coordinates as the nature position orientation
% register all other orientations with the common space coordinate
FLIRT

% register the magnitude from first echo
% apply the transformation to local field map


% calculate the angles of B0 with registered local field maps
% R is transformation matrix from individual image space to common space (FLIRT matrix)
% common space coordinates = R* object image space coordinates

% projections of B0 direction on common space coordinates is:
% R*z_prjs, z_prjs derived from DICOM headers for each orientation
z_prjs_c = R'*z_prjs; % (each orientation has own R and z_prjs)
% R is the rotation matrix from image space to common space

% construct k-space kernel for each orientation
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



%% COSMOS reconstruction with closed-form solution

Phase_tissue = zeros(size(phase_tissue));

for t = 1:size(phase_tissue,4)
    Phase_tissue(:,:,:,t) = fftn(phase_tissue(:,:,:,t));
end

kernel_sum = sum(abs(kernel).^2, 4);

chi_cosmos = real( ifftn( sum(kernel .* Phase_tissue, 4) ./ (eps + kernel_sum) ) ) .* msk;
