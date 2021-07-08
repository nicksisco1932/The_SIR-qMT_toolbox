function [ss_tp, gd_tp, pk_tp, brainmask] = determineTimePoints_NS(vol, brainMask, TR,SaveFlag,VOXflag) 
if length(size(vol)) == 5
    echoFlag = 1;
elseif length(size(vol)) == 4
    echoFlag = 0;
else
    warning('Incorrect Matrix Size! We must be a time series!')
    return;
end


%LCB 7/1/2016
if ~exist('TR','var')
    TR = 1;
    warning('assuming TR is 1 second');
end
%% First mask brain and determine steady state and contrast arrival timepoints
switch echoFlag
    case 1
        [nx,ny,nz,ne,nt] = size(vol);
        I = squeeze(vol(:,:,:,2,:)); %only need first echo dataset for this section
        I_v = reshape(I, [nx*ny*nz nt]);
    otherwise
        [nx,ny,nz,nt] = size(vol);
        ne = 1; %LCB 4/6/2016
        I = vol;
        I_v = reshape(I, [nx*ny*nz nt]);
end


brainmask = brainMask;
%LCB 7/1/2016 % sure I'll leave this but it's not necessary
er_brainmask = zeros(nx,ny,nz);
for z = 1:nz
    er_brainmask(:,:,z) = imerode(imfill(brainmask(:,:,z)),strel('disk', 2, 0)); %AMS changed strel('disk',5,0) to strel('disk',2,0), bigger mask
end
%figure, montage(permute(er_brainmask, [1 2 4 3])); title('Brain mask');
brainmask = er_brainmask;

if isequal(SaveFlag, 1);
    save brainmask.mat brainmask
end

switch VOXflag
    case 0
        %Find steady state and contrast arrival timepoints
        mean_tc = mean(reshape(I, [nx*ny*nz nt])); %mean time course
        slope = floor(diff(mean_tc)); %find slope of curve to first determine steady state location
        if slope(1) < -20 || slope(1) > 20 %include positive slope - important for SAGE
            ss_tp = find(slope >= -1 & slope <= 1, 1, 'first') + round(2/TR); %steady state reached, add 2 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1; %AMS 7/21/2016 -1 <= slope <= 1
            if ss_tp > nt/2
                warning('Did not find steady-state in first half of time-course, trying again')
                ss_tp = find(slope >= -20 & slope <= 20, 1, 'first') + round(2/TR); %steady state reached, add 2 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1; %AMS 7/21/2016 -1 <= slope <= 1
            end
        else
            ss_tp = 1; %assume no dummy scans included in data
        end
        if isempty(ss_tp) %catch empty ss_tp??
            ss_tp = 1;
        end
        [pks,locs,w] = findpeaks(-mean_tc(:,ss_tp:end));
        [~, gd_index] = max(pks);
        gd_tp = locs(gd_index) + ss_tp + 1 - round(3/TR) - floor(w(gd_index)*1.5); %contrast arrival
        pk_tp = locs(gd_index) + ss_tp - round(3/TR) + 2; %time to peak

    case 1
        warning('Needs to be coded. Lucky you!');
        return;
end

end