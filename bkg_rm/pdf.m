function [lfs,mask_ero] = pdf(tfs,mask,voxelSize,ker_rad,weight,theta)


if ~ exist('theta','var') || isempty(theta)
    theta = 0;
end

% add zero slices
tfs = padarray(tfs,[0 0 20]);
mask = padarray(mask,[0 0 20]);
weight = padarray(weight,[0 0 20]);



[nv,np,ns] = size(tfs);


%%%%%%%%%%%%%%%%%%%%% needs scaling with FOV %%%%%%%%%%%%%%%%%%%%%%%
% create K-space filter kernel D

FOV = voxelSize.*[nv,np,ns];
FOVx = FOV(1);
FOVy = FOV(2);
FOVz = FOV(3);

x = -nv/2:nv/2-1;
y = -np/2:np/2-1;
z = -ns/2:ns/2-1;
[kx,ky,kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
D = 1/3 - kz.^2./(kx.^2 + ky.^2 + kz.^2);

D(nv/2+1,np/2+1,ns/2+1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(theta) >= 0.1 % 18/pi = 5.73d
    disp('rotate dipole kernel');
    %% rotate the kernel to match the angle of acqisition
    [X,Y,Z] = ndgrid(x,y,z);
    % rotM = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    rotM = [1,0,0; 0,cos(theta),-sin(theta); 0,sin(theta),cos(theta)];

    afterR = rotM*[X(:)' ; Y(:)' ; Z(:)'];
    finalM = afterR + repmat([nv/2+1;np/2+1;ns/2+1],[1 nv*np*ns]);
    %Index = round(finalM);
    Index = finalM;

    IndexX = Index(1,:);
    IndexY = Index(2,:);
    IndexZ = Index(3,:);

    IndexX(IndexX<=1) = 1;
    IndexY(IndexY<=1) = 1;
    IndexZ(IndexZ<=1) = 1;

    IndexX(IndexX>nv) = nv;
    IndexY(IndexY>np) = np;
    IndexZ(IndexZ>ns) = ns;

    % bilinear interpolation
    D_r = zeros(size(D));
    for i = 1:length(IndexX)
        ix = IndexX(i);
        iy = IndexY(i);
        iz = IndexZ(i);
        d = iz - floor(iz);
        d1 = iy - floor(iy);
        D_r(i) = D(ix,floor(iy),floor(iz))*(1-d)*(1-d1) + D(ix,floor(iy),ceil(iz))*d*(1-d1) + D(ix,ceil(iy),floor(iz))*d1*(1-d) + D(ix,ceil(iy),ceil(iz))*d*d1;
    end

    D = D_r;
end

D = fftshift(D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 1-mask; % region outside the mask
W = weight.*mask; % weights


b = M.*ifftn(D.*fftn(W.*W.*tfs));
b = b(:);

res = cgs(@Afun,b,1e-6, 200);
lfs = mask.*real(tfs-ifftn(D.*fftn(M.*reshape(res,[nv,np,ns]))));

function y = Afun(x)
    x = reshape(x,[nv,np,ns]);
    y = M.*ifftn(D.*fftn(W.*W.*ifftn(D.*fftn(M.*x))));
    y = y(:);
end


% erode the edge
rx = round(ker_rad/voxelSize(1));
ry = round(ker_rad/voxelSize(2));
rz = round(ker_rad/voxelSize(3));
[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 < 1);
ker = h/sum(h(:));
csh = [rx,ry,rz];
imsize = size(mask);
mask_tmp = circshift(real(ifftn(fftn(mask).*fftn(ker,imsize))),-csh);
mask_ero = zeros(imsize);
mask_ero(mask_tmp > 1-1/sum(h(:))) = 1; % no error tolerance
lfs = lfs.*mask_ero;


% remove added zero slices
lfs = lfs(:,:,21:end-20);
mask_ero = mask_ero(:,:,21:end-20);

end
