function DicomWrite(FileName,Image,VS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A function to create a tiff stack and a Dicom File containing the images from a matlab 3D matrix
%Images indexed by 3rd dimension show 1st and 2nd dimension
%
%Pass to the function: 
%                     Image (N1xN2xN3 matrix)
%                     Filename (with or without .tif extension - either
%                     will work
%                     VS - Voxel Size in mm either scalar (for isotropic) or
%                     3-vector. If not specified, defaults to 1 mm
%                     isotropic
%Peter Niedbalski, 2018
%Functioning as of 4/2/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3
    VS = [1 1 1];
end

if numel(VS)==1
   VS = [VS(1) VS(1) VS(1)]; 
elseif numel(VS) ~= 3
   error('VS needs to be either a scalar or 3-d vector.'); 
end

%Normalize Image
Image = abs(Image)/max(max(max(abs(Image))));
%I was running into issues with the format that I saved with. uint16 seems
%to work well, so I need to convert my double image to uint16
Image = uint16(round(abs(Image)*65535));

FileNamedcm = 'DicomImage';

%Make sure that filename has .tif extension
if strcmp(FileName((length(FileName)-3):length(FileName)),'.tif')~=1
    FileNamedcm = [FileName '.dcm'];
    FileName = [FileName '.tif'];
else
    FileNamedcm = [FileName(1:(length(FileName)-4)) '.dcm'];
end

% create a dicom header with the relevant information (Copied from
% dicomwritevolume function from online)
info.SliceThickness = VS(3);
info.SpacingBetweenSlices = VS(3); %ITK-Snap fix for 3D images 23Jul2019-ba
info.ImagerPixelSpacing = VS(1:2);
info.PixelSpacing = VS(1:2);
info.Width = size(Image,1);
info.Height = size(Image,2);
info.ColorType = 'grayscale';
info.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.4'; % MR image storage
info.TransferSyntaxUID = '1.2.840.10008.1.2.1'; % Explicit VR Little Endian
info.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4'; % MR Image Storage
info.PhotometricInterpretation = 'MONOCHROME2';
info.PixelRepresentation = 0;
info.WindowCenter = 25000; %These are the parameters to change to fix the window and level of the dicom image
info.WindowWidth = 50000;
info.RescaleIntercept = -1024;
info.RescaleSlope = 1;
info.RescaleType = 'HU';
%Write Dicom File
dicomwrite(reshape(Image,size(Image,1),size(Image,2),1,size(Image,3)),FileNamedcm,info,'CreateMode','copy');

%Create Tiff File and write image
% t = Tiff(FileName,'w');
% tagstruct.ImageLength = size(Image,1);
% tagstruct.ImageWidth = size(Image,2);
% tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
% tagstruct.BitsPerSample = 16;
% tagstruct.SamplesPerPixel = 1;
% tagstruct.RowsPerStrip = 1;
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% tagstruct.Software = 'MATLAB';
% tagstruct.SampleFormat = 1;
% %tagstruct.SubIFD = size(Image,3)-1 ;
% setTag(t,tagstruct)
% 
% %Loop through image and save images in subdirectories
% for i = 1:size(Image,3)
%     write(t,squeeze(Image(:,:,i)),-1);
%     if i < size(Image,3)
%         writeDirectory(t);
%         tagstruct.ImageLength = size(Image,1);
%         tagstruct.ImageWidth = size(Image,2);
%         tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%         tagstruct.BitsPerSample = 16;
%         tagstruct.SamplesPerPixel = 1;
%         tagstruct.RowsPerStrip = 1;
%         tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%         tagstruct.Software = 'MATLAB';
%         tagstruct.SampleFormat = 1;
%         setTag(t,tagstruct)
%     end
% end
% close(t)