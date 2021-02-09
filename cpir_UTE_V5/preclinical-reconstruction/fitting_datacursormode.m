function fitting_datacursormode(RawDataX,RawDataY,FitVal,FitFunc,FitOpts)

myFig = imslice(FitVal,'Fitting Result');

NCmap = jet;
NCmap(1,:) = [0 0 0];
colormap(myFig,NCmap)

rawfig = figure('Name','Data from which Fit is Derived');

dcm_obj = datacursormode(myFig);
set(dcm_obj,'UpdateFcn',{@fittingdisplayfunc,myFig,rawfig,RawDataX,RawDataY,FitFunc,FitOpts})


