function FID = splitFID(FIDs,NTE)
%pass FIDs, trajectories, and number of TE times to separate FID and
%trajectories into 3-D and 4-D matrices containing each individual image
%FID and Trajectories

npts = size(FIDs,1);
TotRays = size(FIDs,2);
RaysperImage = TotRays/NTE;

FID = zeros(npts,RaysperImage,NTE);
%Traj = zeros(3,npts,RaysperImage,NTE);

RaysofImage = 1:NTE:TotRays;

for i = 1:NTE
    FID(:,:,i) = FIDs(:,RaysofImage+(i-1));
    %Traj(:,:,:,i) = traj(:,:,RaysofImage+(i-1));
end
