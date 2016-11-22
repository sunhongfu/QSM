% run QSM individually on each orientation




% define common space coordinates as the nature position orientation
% register all other orientations with the common space coordinate
unix('/usr/share/fsl/5.0/bin/flirt -in /home/hongfu/Desktop/9175/QSM_SPGRE_LEFT/QSM_SPGR_GE/src/mag1.nii -ref /home/hongfu/Desktop/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/src/mag1.nii -out /home/hongfu/Desktop/9175/QSM_SPGRE_LEFT/mag_flirt.nii -omat /home/hongfu/Desktop/9175/QSM_SPGRE_LEFT/mag_flirt.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear');

unix('/usr/share/fsl/5.0/bin/flirt -in /home/hongfu/Desktop/9175/QSM_SPGRE_RIGHT/QSM_SPGR_GE/src/mag1.nii -ref /home/hongfu/Desktop/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/src/mag1.nii -out /home/hongfu/Desktop/9175/QSM_SPGRE_RIGHT/mag_flirt.nii -omat /home/hongfu/Desktop/9175/QSM_SPGRE_RIGHT/mag_flirt.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear');

unix('/usr/share/fsl/5.0/bin/flirt -in /home/hongfu/Desktop/9175/QSM_SPGRE_FORWARD/QSM_SPGR_GE/src/mag1.nii -ref /home/hongfu/Desktop/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/src/mag1.nii -out /home/hongfu/Desktop/9175/QSM_SPGRE_FORWARD/mag_flirt.nii -omat /home/hongfu/Desktop/9175/QSM_SPGRE_FORWARD/mag_flirt.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear');

unix('/usr/share/fsl/5.0/bin/flirt -in /home/hongfu/Desktop/9175/QSM_SPGRE_BACKWARD/QSM_SPGR_GE/src/mag1.nii -ref /home/hongfu/Desktop/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/src/mag1.nii -out /home/hongfu/Desktop/9175/QSM_SPGRE_BACKWARD/mag_flirt.nii -omat /home/hongfu/Desktop/9175/QSM_SPGRE_BACKWARD/mag_flirt.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear');

% apply the transformation to local field map
unix('/usr/share/fsl/5.0/bin/flirt -in /home/hongfu/Desktop/9175/QSM_SPGRE_LEFT/QSM_SPGR_GE/RESHARP/lfs_resharp_tik_0.0001_num_500.nii -applyxfm -init /home/hongfu/Desktop/9175/QSM_SPGRE_LEFT/mag_flirt.mat -out /home/hongfu/Desktop/9175/QSM_SPGRE_LEFT/QSM_SPGR_GE/RESHARP/lfs_resharp_flirt.nii -paddingsize 0.0 -interp trilinear -ref /home/hongfu/Desktop/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/RESHARP/lfs_resharp_tik_0.0001_num_500.nii');

unix('/usr/share/fsl/5.0/bin/flirt -in /home/hongfu/Desktop/9175/QSM_SPGRE_RIGHT/QSM_SPGR_GE/RESHARP/lfs_resharp_tik_0.0001_num_500.nii -applyxfm -init /home/hongfu/Desktop/9175/QSM_SPGRE_RIGHT/mag_flirt.mat -out /home/hongfu/Desktop/9175/QSM_SPGRE_RIGHT/QSM_SPGR_GE/RESHARP/lfs_resharp_flirt.nii -paddingsize 0.0 -interp trilinear -ref /home/hongfu/Desktop/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/RESHARP/lfs_resharp_tik_0.0001_num_500.nii');

unix('/usr/share/fsl/5.0/bin/flirt -in /home/hongfu/Desktop/9175/QSM_SPGRE_FORWARD/QSM_SPGR_GE/RESHARP/lfs_resharp_tik_0.0001_num_500.nii -applyxfm -init /home/hongfu/Desktop/9175/QSM_SPGRE_FORWARD/mag_flirt.mat -out /home/hongfu/Desktop/9175/QSM_SPGRE_FORWARD/QSM_SPGR_GE/RESHARP/lfs_resharp_flirt.nii -paddingsize 0.0 -interp trilinear -ref /home/hongfu/Desktop/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/RESHARP/lfs_resharp_tik_0.0001_num_500.nii');

unix('/usr/share/fsl/5.0/bin/flirt -in /home/hongfu/Desktop/9175/QSM_SPGRE_BACKWARD/QSM_SPGR_GE/RESHARP/lfs_resharp_tik_0.0001_num_500.nii -applyxfm -init /home/hongfu/Desktop/9175/QSM_SPGRE_BACKWARD/mag_flirt.mat -out /home/hongfu/Desktop/9175/QSM_SPGRE_BACKWARD/QSM_SPGR_GE/RESHARP/lfs_resharp_flirt.nii -paddingsize 0.0 -interp trilinear -ref /home/hongfu/Desktop/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/RESHARP/lfs_resharp_tik_0.0001_num_500.nii');

% calculate the angles of B0 with registered local field maps
% R is transformation matrix from individual image space to common space (FLIRT matrix)
% common space coordinates = R* object image space coordinates
load QSM_SPGRE_LEFT/mag_flirt.mat -ASCII
R_t(:,:,1) = mag_flirt(1:3,1:3);
load QSM_SPGRE_RIGHT/mag_flirt.mat -ASCII
R_t(:,:,2) = mag_flirt(1:3,1:3);
load QSM_SPGRE_FORWARD/mag_flirt.mat -ASCII
R_t(:,:,3) = mag_flirt(1:3,1:3);
load QSM_SPGRE_BACKWARD/mag_flirt.mat -ASCII
R_t(:,:,4) = mag_flirt(1:3,1:3);
R_t(:,:,5) = eye(3);

% (each orientation has own R and z_prjs)
% R is the rotation matrix from image space to common space
load QSM_SPGRE_LEFT/QSM_SPGR_GE/z_prjs.mat
z_prjs_o(:,1) = z_prjs';
load QSM_SPGRE_RIGHT/QSM_SPGR_GE/z_prjs.mat
z_prjs_o(:,2) = z_prjs';
load QSM_SPGRE_FORWARD/QSM_SPGR_GE/z_prjs.mat
z_prjs_o(:,3) = z_prjs';
load QSM_SPGRE_BACKWARD/QSM_SPGR_GE/z_prjs.mat
z_prjs_o(:,4) = z_prjs';
load QSM_SPGRE_CENTER/QSM_SPGR_GE/z_prjs.mat
z_prjs_o(:,5) = z_prjs';

for i = 1:5
    z_prjs_c(:,i) = R_t(:,:,i)'*z_prjs_o(:,i);
end


%% COSMOS reconstruction with closed-form solution
% load in registered local field shift maps
unix('gunzip -f /home/hongfu/Desktop/9175/QSM_SPGRE_LEFT/QSM_SPGR_GE/RESHARP/*.gz');
nii = load_nii('/home/hongfu/Desktop/9175/QSM_SPGRE_LEFT/QSM_SPGR_GE/RESHARP/lfs_resharp_flirt.nii');
lfs(:,:,:,1) = double(nii.img);
unix('gunzip -f /home/hongfu/Desktop/9175/QSM_SPGRE_RIGHT/QSM_SPGR_GE/RESHARP/*.gz');
nii = load_nii('/home/hongfu/Desktop/9175/QSM_SPGRE_RIGHT/QSM_SPGR_GE/RESHARP/lfs_resharp_flirt.nii');
lfs(:,:,:,2) = double(nii.img);
unix('gunzip -f /home/hongfu/Desktop/9175/QSM_SPGRE_FORWARD/QSM_SPGR_GE/RESHARP/*.gz');
nii = load_nii('/home/hongfu/Desktop/9175/QSM_SPGRE_FORWARD/QSM_SPGR_GE/RESHARP/lfs_resharp_flirt.nii');
lfs(:,:,:,3) = double(nii.img);
unix('gunzip -f /home/hongfu/Desktop/9175/QSM_SPGRE_BACKWARD/QSM_SPGR_GE/RESHARP/*.gz');
nii = load_nii('/home/hongfu/Desktop/9175/QSM_SPGRE_BACKWARD/QSM_SPGR_GE/RESHARP/lfs_resharp_flirt.nii');
lfs(:,:,:,4) = double(nii.img);
nii = load_nii('/home/hongfu/Desktop/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/RESHARP/lfs_resharp_tik_0.0001_num_500.nii');
lfs(:,:,:,5) = double(nii.img);

mask = and(and(and(and(lfs(:,:,:,1),lfs(:,:,:,2)),lfs(:,:,:,3)),lfs(:,:,:,4)),lfs(:,:,:,5));
mask = double(mask);

% construct k-space kernel for each orientation
% create K-space filter kernel D
%%%%% make this a seperate function in the future
vox = [1 1 1];

Nx = size(lfs,1);
Ny = size(lfs,2);
Nz = size(lfs,3);

for i = 1:5
FOV = vox.*[Nx,Ny,Nz];
FOVx = FOV(1);
FOVy = FOV(2);
FOVz = FOV(3);

x = -Nx/2:Nx/2-1;
y = -Ny/2:Ny/2-1;
z = -Nz/2:Nz/2-1;
[kx,ky,kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
% D = 1/3 - kz.^2./(kx.^2 + ky.^2 + kz.^2);
D = 1/3 - (kx.*z_prjs_c(1,i)+ky.*z_prjs_c(2,i)+kz.*z_prjs_c(3,i)).^2./(kx.^2 + ky.^2 + kz.^2);
D(floor(Nx/2+1),floor(Ny/2+1),floor(Nz/2+1)) = 0;
D = fftshift(D);

kernel(:,:,:,i) = D;
end


for i = 1:5
    lfs_k(:,:,:,i) = fftn(lfs(:,:,:,i).*mask);
end

kernel_sum = sum(abs(kernel).^2, 4);

chi_cosmos = real( ifftn( sum(kernel .* lfs_k, 4) ./ (eps + kernel_sum) ) ) .* mask;

nii = make_nii(chi_cosmos,vox);
save_nii(nii,'cosmos_5.nii');


%%%% inversion of each orientation on registered LFS
nii = load_nii('/media/data/QSM_data/COSMOS_renzo/sorted/patient/9175/QSM_SPGRE_CENTER/QSM_SPGR_GE/src/mag1.nii');
mag = double(nii.img);
for i = 1:5
    disp('--> TV susceptibility inversion on RESHARP...');
    sus_resharp(:,:,:,i) = tvdi(lfs(:,:,:,i),mask,vox,tv_reg,mag,z_prjs_c(:,i),inv_num); 
   
    % save nifti
    nii = make_nii(sus_resharp(:,:,:,i).*mask,vox);
    save_nii(nii,['sus_resharp_tik_', num2str(tik_reg), '_tv_', num2str(tv_reg), '_num_', num2str(inv_num), '_orien_' num2str(i) '.nii']);
end
