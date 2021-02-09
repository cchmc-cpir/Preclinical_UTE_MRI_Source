function [Image_Out,kspace,Recon_Params] = base_noncart_recon2D(ImageSize,data,traj)
%% A Function written to reconstruct Images when K-space data and trajectories are passed to it
% Uses Scott Robertson's reconstruction code - This just makes it more
% modular and easy to implement - This is for 3D data
% 
% ImageSize - Output image matrix size -If Scalar, Isotropic, ImageSize
%                                       If 3-Vector,
%
% data - KSpace Data in column vector (N x 1)
%
% traj - point in kspace corresponding to the data vector - columns for
% x,y, and z. (N x 3)

%Make sure ImageSize is Meaningful
if numel(ImageSize)==1
   ImageSize = [ImageSize(1) ImageSize(1) ImageSize(1)]; 
elseif numel(ImageSize) ~= 3
   error('ImageSize needs to be either a scalar or 3-d vector.'); 
end

kernel.sharpness = .3;
kernel.extent = 6*kernel.sharpness;
overgrid_factor = 3;
output_image_size = ImageSize;
nDcfIter = 10;
deapodizeImage = true; %set to true -- BA_06Aug2019
nThreads = 12;
cropOvergriddedImage = true();
verbose = true();

%% Save the important parameters
Recon_Params.KernelSharpness = kernel.sharpness;
Recon_Params.KernelExtent = kernel.extent;
Recon_Params.OverGridding = overgrid_factor;
Recon_Params.DCIterations = nDcfIter;
Recon_Params.nThreads = nThreads;
Recon_Params.Recon_Func = 'base_noncart_recon.m';
Recon_Params.Deapodize = deapodizeImage;

%  Choose kernel, proximity object, and then create system model
kernelObj = Recon.SysModel.Kernel.Gaussian(kernel.sharpness, kernel.extent, verbose);
%kernelObj = Recon.SysModel.Kernel.KaiserBessel(kernel.sharpness, kernel.extent, verbose);
%kernelObj = Recon.SysModel.Kernel.Sinc(kernel.sharpness, kernel.extent, verbose);

proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);
%proxObj = Recon.SysModel.Proximity.L1Proximity(kernelObj, verbose);
clear kernelObj;
systemObj = Recon.SysModel.MatrixSystemModel(traj, overgrid_factor, ...
    output_image_size, proxObj, verbose);

% Choose density compensation function (DCF)
dcfObj = Recon.DCF.Iterative(systemObj, nDcfIter, verbose);
%dcfObj = Recon.DCF.Voronoi(traj, header, verbose);
%dcfObj = Recon.DCF.Analytical3dRadial(traj, verbose);
%dcfObj = Recon.DCF.Unity(traj, verbose);

% Choose Reconstruction Model
reconObj = Recon.ReconModel.LSQGridded(systemObj, dcfObj, verbose);
clear modelObj;
clear dcfObj;
reconObj.crop = cropOvergriddedImage;
reconObj.deapodize = deapodizeImage;
% reconObj.apodizeKspace = apodizeKspace %BA 06Aug2019

% Reconstruct image using trajectories in pixel units
[Image_Out,kspace] = reconObj.reconstruct(data, traj); % BA 06Aug2019 add in section that runs apodization on kspace; add in a read in variable
Image_Out = rot90(flip(Image_Out));
% for i = 1:ImageSize
%     Image_Out(:,:,i) = fliplr(rot90(Image_Out(:,:,i),1));
% end
% %imslice(abs(Image_Out))
% for i = 1:ImageSize
%     test = reshape(Image_Out(i,:,:),ImageSize(2),ImageSize(3));
%     test = fliplr(test);
%     Image_Out(i,:,:) = reshape(test,1,ImageSize,ImageSize);
% end
