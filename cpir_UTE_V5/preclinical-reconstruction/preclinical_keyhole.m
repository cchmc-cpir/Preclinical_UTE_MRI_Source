function KeyImages = preclinical_keyhole(fid,traj,Method_Params)
    ImSize = Method_Params.MatrixSize;
    %Get number pof projections per trigger
    PPT = Method_Params.ProjPerTrig;
    KeyImages = zeros(ImSize(1),ImSize(2),ImSize(3),PPT);
    %Get the keyhole radius - make it a minimum of 1/4 the total k-space radius so that we have adequate sampling 
    KeyRad = max(floor(size(fid,1)/4),floor(sqrt(size(fid,2)/4/pi/PPT)));
    Keyhole = fid;
    Keyhole(1:KeyRad,:) = NaN;
    KeyKSpace = zeros(size(fid,1),size(fid,2),PPT);
    for ii = 1:PPT
        KeyKSpace(:,:,ii) = Keyhole;
        keyk0 = abs(fid(1,ii:PPT:size(fid,2)));
        MeanKeyk0 = mean(keyk0);
        for jj = 1:size(fid,2)
            KeyKSpace(:,jj,ii) = KeyKSpace(:,jj,ii)*MeanKeyk0/abs(fid(1,jj));
        end
        KeyKSpace(1:KeyRad,ii:PPT:size(fid,2),ii) = fid(1:KeyRad,ii:PPT:size(fid,2));
        keyfid = reshape(squeeze(KeyKSpace(:,:,ii)),1,[])';
        findzeros = find(isnan(keyfid));
        ThisKeyTraj = traj;
        ThisKeyTraj(findzeros,:) = [];
        keyfid(findzeros) = [];

        %Remove Points that are beyond a trajectory radius of 0.5
        traj_recon_1 = ThisKeyTraj;
        rad = sqrt(traj_recon_1(:,1).^2+traj_recon_1(:,2).^2+traj_recon_1(:,3).^2);
        del_pts = find(rad > 0.5);
        keyfid(del_pts) = [];
        traj_recon_1(del_pts,:) = [];

        [KeyImages(:,:,:,ii),~,~] = base_noncart_recon(ImSize,keyfid,traj_recon_1);
    end