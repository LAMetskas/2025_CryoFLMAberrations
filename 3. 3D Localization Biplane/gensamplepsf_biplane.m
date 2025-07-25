function [samplepsf,startx,starty,startz, dz,dx] = gensamplepsf_biplane(psfobj,pixelsize,psfsize,boxsize,bin,Nzs, dist_biplane)
zpos = linspace(-1.3,1.3,Nzs)';
xpos = zeros(Nzs,1);
ypos = zeros(Nzs,1);
I = ones(Nzs,2);
bg = zeros(Nzs,2);
x1 = cat(2,xpos,ypos,zpos,I,bg);
w = [1,1,1,1,1];


psfobj.Pixelsize = pixelsize/bin;
psfobj.PSFsize = psfsize;
psfobj.Boxsize = boxsize*bin;

[psf2d_fit] = genpsf_biplane_real(x1,w,psfobj,dist_biplane);


psf_fit = [];
for ii = 1:2
    psf_fit = cat(2,psf_fit,squeeze(psf2d_fit(:,:,:,ii)));
end
samplepsf = cell(2,1);
for ii = 1:2
    f0 = psf2d_fit(:,:,:,ii);
    N = size(f0,1);
    Nz = size(f0,3);
    F = reshape(permute(f0,[2,1,3]),N*N*Nz,1);

    samplepsf{ii} = F;
end
startx = -0.5*psfobj.Pixelsize*psfobj.Boxsize;
starty = -0.5*psfobj.Pixelsize*psfobj.Boxsize;
startz = zpos(1);
dz = zpos(2)-zpos(1);
dx = psfobj.Pixelsize;
end
