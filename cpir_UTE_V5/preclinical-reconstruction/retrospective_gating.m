function Gating_Vector = retrospective_gating(FID,Method_Params)

% Loop through different echoes
Gating_Vector = zeros(1,size(FID,2),size(FID,3));

for i = 1:size(FID,3)
    %Isolate k0
    k0 = abs(FID(1,:,i));
    %Normalize k0 (since it's easier)
    k0 = k0/max(k0);
    %smooth k0 (Was 25)
    sm_k0 = smooth(k0,25);
    %Take first derivative
    D1 = diff(sm_k0,1);
    %Plug in the last point so our arrays remain the same length
    D1(length(k0)) = D1(length(k0)-1);
    %Smooth first derivative
    sm_D1 = smooth(D1,10);
    %Take second derivative
    D2 = diff(sm_D1);
    %Plug in the last point so our arrays remain the same length
    D2(length(k0)) = D2(length(k0)-1);
%% Select projections at End Expiration
    explim = 0.01 * max(sm_D1);
    
    Gate_Vec = zeros(1,length(k0));
    if sm_D1(1) > -explim && sm_D1(1) < explim
        Gate_Vec(1) = 1;
    end
    last_val = 1;
    for j = 3:size(FID,2)
        if sm_D1(j) > -explim && sm_D1(j) < explim
            if Gate_Vec(j-1) == 0
                Gate_Vec(j) = -last_val;
            else
                Gate_Vec(j) = Gate_Vec(j-1);
            end
            last_val = Gate_Vec(j);
        end
    end     
    allvec = find(Gate_Vec ~= 0);
    
    med = (max(sm_k0) - min(sm_k0))/2 + min(sm_k0);
    
    for j = allvec
        if sm_k0(j) > med
            Gate_Vec(j) = 1;
        elseif sm_k0(j) < med
            Gate_Vec(j) = -1;
        end
    end
    
    Gate_1 = find(Gate_Vec == 1);
    Gate_2 = find(Gate_Vec == -1);
    
    if length(Gate_1) > length(Gate_2)
        Exp_Vec = Gate_1;
    else
        Exp_Vec = Gate_2;
    end
%% Select Projections for End Inspiration
    insplim = 0.05 * max(sm_D1);
    
    Gate_Vec = zeros(1,length(k0));
    if sm_D1(1) > -insplim && sm_D1(1) < insplim
        Gate_Vec(1) = 1;
    end
    last_val = 1;
    for j = 3:size(FID,2)
        if sm_D1(j) > -insplim && sm_D1(j) < insplim
            if Gate_Vec(j-1) == 0
                Gate_Vec(j) = -last_val;
            else
                Gate_Vec(j) = Gate_Vec(j-1);
            end
            last_val = Gate_Vec(j);
        end
    end     
    allvec = find(Gate_Vec ~= 0);
    
    med = (max(sm_k0) - min(sm_k0))/2 + min(sm_k0);
    
    for j = allvec
        if sm_k0(j) > med
            Gate_Vec(j) = 1;
        elseif sm_k0(j) < med
            Gate_Vec(j) = -1;
        end
    end
    
    Gate_1 = find(Gate_Vec == 1);
    Gate_2 = find(Gate_Vec == -1);
    
    if length(Gate_1) > length(Gate_2)
        Insp_Vec = Gate_2;
    else
        Insp_Vec = Gate_1;
    end
    
    Gating_Vector(1,Exp_Vec,i) = 1;
    Gating_Vector(1,Insp_Vec,i) = -1;
    
    %% Display Retro Gating Results
    figure('Name',['Gating for Image ' num2str(i)])
    Gate_1 = find(Gating_Vector(:,:,i) == 1);
    Gate_2 = find(Gating_Vector(:,:,i) == -1);
    Index = 1:length(k0);
    time_axis = (Index-1)*Method_Params.TR;
    plot(time_axis,sm_k0,'ko')
    hold on
    plot(time_axis(Gate_1),sm_k0(Gate_1),'ro','MarkerFaceColor','r')
    plot(time_axis(Gate_2),sm_k0(Gate_2),'bo','MarkerFaceColor','b')
    hold off
    legend('All Data','End Expiration','End Inspiration')
end



