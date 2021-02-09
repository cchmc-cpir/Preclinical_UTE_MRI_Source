function Image = washout_recon(fid,traj,Method_Params,Config_Params)

%Function to reconstruct washout images acquired on the Bruker

%Let's do things all keyhole-like
NIm = Config_Params.NWashin + Config_Params.NWashout;

%I'll want to be more clever about this eventually, but to try to get
%things to work, start with hardcoding some things
fid = fid(:,187:end);

%Figure out what my keyhole radius should be
NPro = size(fid,2);

PPT = Method_Params.ProjPerTrig;

%Find the number of projections between images of the same type
nbtw = PPT*NIm;

%Get the number of projections per image
nperim = round(NPro/NIm);

WI_ind = 1:(Config_Params.NWashin*PPT);
WO_ind = (Config_Params.NWashin*PPT+1):(Config_Params.NWashin*PPT+Config_Params.NWashout*PPT);

count = 1;
while max(WI_ind)<NPro || max(WO_ind)<NPro
    WI_ind = [WI_ind ((Config_Params.NWashin*PPT+Config_Params.NWashout*PPT)*count+1):((Config_Params.NWashin*PPT+Config_Params.NWashout*PPT)*count+Config_Params.NWashin*PPT)];
    WO_ind = [WO_ind ((Config_Params.NWashin*PPT+Config_Params.NWashout*PPT)*count+Config_Params.NWashin*PPT+1):((Config_Params.NWashin*PPT+Config_Params.NWashout*PPT)*(count+1))];
    count = count +1;
end
    
WI_ind(WI_ind>NPro) = [];
WO_ind(WO_ind>NPro) = [];
%Get keyhole radius:
keyrad = round(sqrt(nperim/4/pi));
%keyrad = size(fid,1);

%I think ultimately, we'll want to stack keyhole on keyhole to get RF decay
%and T1 and Washout and whatnot, but start with the simple method first

%Empty out Keyhole

fid_in = fid(:,WI_ind);
keyhole_in = fid(:,WI_ind);
keyhole_in(1:keyrad,:) = nan;

figure('Name','Wash In FIDs')
imagesc(abs(fid_in))
colormap(gray)

fid_out = fid(:,WO_ind);
keyhole_out = fid(:,WO_ind);
keyhole_out(1:keyrad,:) = nan;

figure('Name','Wash Out FIDs')
imagesc(abs(fid_out))
colormap(gray)

traj_in = traj(:,:,WI_ind);
traj_out = traj(:,:,WO_ind);

trajx = reshape(traj_in(1,:,:),1,[])';
trajy = reshape(traj_in(2,:,:),1,[])';
trajz = reshape(traj_in(3,:,:),1,[])';

trajC_in = [trajx trajy trajz];

trajx = reshape(traj_out(1,:,:),1,[])';
trajy = reshape(traj_out(2,:,:),1,[])';
trajz = reshape(traj_out(3,:,:),1,[])';

trajC_out = [trajx trajy trajz];

indices = 1:PPT;
nbtw_in = PPT*Config_Params.NWashin;
indx_in = [];
in_trigs = size(fid_in,2)/PPT;
for i = 1:in_trigs
    indx_in = [indx_in (i-1)*nbtw_in+indices];
end

indices = 1:PPT;
nbtw_out = PPT*Config_Params.NWashout;
indx_out = [];
out_trigs = size(fid_out,2)/PPT;
for i = 1:out_trigs
    indx_out = [indx_out (i-1)*nbtw_out+indices];
end

trajx = reshape(traj(1,:,:),1,[])';
trajy = reshape(traj(2,:,:),1,[])';
trajz = reshape(traj(3,:,:),1,[])';

trajC = [trajx trajy trajz];

ImSize = Method_Params.MatrixSize;
Image = zeros(ImSize(1),ImSize(2),ImSize(3),NIm); 

for i = 1:Config_Params.NWashin
    key = keyhole_in;
    keyind = indx_in+PPT*(i-1);
    keyind(keyind>size(keyhole_in,2)) = [];
    key(1:keyrad,keyind) = fid_in(1:keyrad,keyind);
    
%     figure('Name','FID Check 1')
%     imagesc(abs(key))
%     colormap(gray)
    
    meank0 = mean(abs(key(1,keyind)));
    for j = 1:size(key,2)
        key(:,j) = key(:,j)*meank0/fid_in(1,j);
    end
    
%     figure('Name','FID Check 2')
%     imagesc(abs(key))
%     colormap(gray)
    
    keyC = reshape(key,1,[])';
    
    keytraj = trajC_in;
    
    nanind = isnan(keyC);
    keyC(nanind) = [];
    keytraj(nanind,:) = [];
    
    rad = sqrt(keytraj(:,1).^2+keytraj(:,2).^2+keytraj(:,3).^2);
    del_pts = find(rad > 0.5);
    
    keyC(del_pts) = [];
    keytraj(del_pts,:) = [];
    
    [Image(:,:,:,i)] = base_noncart_recon(ImSize,keyC,keytraj);
end

for i = 1:Config_Params.NWashout
    key = keyhole_out;
    keyind = indx_out+PPT*(i-1);
    keyind(keyind>size(keyhole_out,2)) = [];
    key(1:keyrad,keyind) = fid_out(1:keyrad,keyind);
    
    figure('Name','FID Check 1')
    imagesc(abs(key))
    colormap(gray)
    
    meank0 = mean(abs(key(1,keyind)));
    for j = 1:size(key,2)
        key(:,j) = key(:,j)*meank0/fid_out(1,j);
    end
    
    figure('Name','FID Check 2')
    imagesc(abs(key))
    colormap(gray)
    
    keyC = reshape(key,1,[])';
    
    keytraj = trajC_out;
    
    nanind = isnan(keyC);
    keyC(nanind) = [];
    keytraj(nanind,:) = [];
    
    rad = sqrt(keytraj(:,1).^2+keytraj(:,2).^2+keytraj(:,3).^2);
    del_pts = find(rad > 0.5);
    
    keyC(del_pts) = [];
    keytraj(del_pts,:) = [];
    
    [Image(:,:,:,i+Config_Params.NWashin)] = base_noncart_recon(ImSize,keyC,keytraj);
end
test = 1;

    

