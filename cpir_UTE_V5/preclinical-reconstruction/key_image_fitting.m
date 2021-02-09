function [M0,FA_Map] = key_image_fitting(Key_Images,Mask)
disp('Beginning Key Fitting')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to generate flip angle map based on Key images reconstructed from
%HP Xenon (in vivo) data
%Going to assume T1 is infinite (since on the timescale of a breath, it
%basically is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Im1 = abs(squeeze(Key_Images(:,:,:,1)));

Key_Images = abs(Key_Images)/max(abs(Key_Images(:)));

meanval = mean(Im1(Mask==1));

%Generate Fitting Function
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0 0 0],...
               'Upper',[Inf 90 Inf],...
               'StartPoint',[meanval 27 .01]);
fiteq =fittype('a.*sind(b).*(cosd(b).^(x-1))+c','options',fo);

xdata = (1:size(Key_Images,4))';

numFits = nnz(Mask);

counter = 1;
%Preallocate memory
M0 = zeros(size(Key_Images,1),size(Key_Images,2),size(Key_Images,3));
FA_Map = zeros(size(Key_Images,1),size(Key_Images,2),size(Key_Images,3));
for ii = 1:size(Key_Images,1)
    for jj = 1:size(Key_Images,2)
        for kk = 1:size(Key_Images,3)
            if Mask(ii,jj,kk) == 1
                ydata = abs(squeeze(Key_Images(ii,jj,kk,:)));
                fo.StartPoint(1) = ydata(1);
                KeyFit = fit(xdata,ydata,fiteq,fo);
                M0(ii,jj,kk) = KeyFit.a;
                FA_Map(ii,jj,kk) = KeyFit.b;
                if counter == numFits/4
                    disp('25% Done with Key Fitting')
                elseif counter == numFits/2
                    disp('50% Done with Key Fitting')
                elseif counter == numFits/4*3
                    disp('75% Done with Key Fitting')
                end
                counter = counter + 1;
            end
        end
    end
end
%Display Flip Angle Map in my Fitting Viewer Function
fitting_datacursormode(xdata,abs(Key_Images),FA_Map,fiteq,fo);