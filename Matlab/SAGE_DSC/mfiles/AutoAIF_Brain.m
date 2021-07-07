%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            Automatic AIF Determination                                                   %
% This function determines the delta R2 star curve for the AIF - note that it does not convert to concentration.           %
% Determining AIF using delta R1 is currently not robust. 12/20/2015                                                       %
%                                                                                                                          %
%                                                                                                                          %
%                                                                                                                          %
% Code from Ashley Stokes and Jack (?)                                                                                     %
% Edited by Laura Bell 12/18/2015                                                                                          %
%                                                                                                                          %
% Usage: AIF = AutoAIF_Brain(vol, paramsTime, flip, DCEflag);                                                              %
% Inputs:                                                                                                                  %
% - vol is a 4D [nx ny nz nt] or 5D [nx ny nz ne nt] matrix                                                                %
% - paramsTime is an array with TR followed by all TEs in seconds                                                          %
% - flip is in degrees (if only single echo data this variable isn't used - put anything)                                  %
% - DCEflag = 0; no R1 calculations                                                                                        %
%           = 1; used for dual echo data without T1 mapping, using Shmainda's et al. method (patent US 8,670,602 B2)       %
%           = 2; multiecho T1 mapping - currently not coded                                                                %
% Outputs:                                                                                                                 %
% - output will be the AIF array in delta R2star and delta R1 (if DCEflag set) -- NOT concentration!                       %
% - saved outputs in working directory:                                                                                    %
%       brainmask.mat is a binary mask of the brain using SI threshold and region growing                                  %
%       AIFmask.mat is a binary mask of the determined AIF pixels                                                          %
%       AutoAIF_steps.fig is a file with 3 figures:                                                                        %
%           1) signal intentisy threshold to remove background signal                                                      %
%           2) Montage of brain mask slices to be used                                                                     %
%           3) Line graph of the mean SI intensity of masked brain with time points detected                               %
%       AIF_dR2s_timecourse.fig                                                                                            %
%       AIFlocations.fig is a montage with brain slices at peak contrast arrival and AIF locations                         %
%                                                                                                                          %
% Additional functions needed to run this function:                                                                        %
% - Region Growing by Christian Wuerslin                                                                                   %
% http://www.mathworks.com/matlabcentral/fileexchange/41666-fast-3d-2d-region-growing--mex-                                %
%                                                                                                                          %
% Edits:                                                                                                                   %
% - AMS 2/2/2016                                                                                                           %
%       Fixed taking the top 10 pixels if there are less than 10 pixels for the AIF                                        %
% - AMS 6/20/2016                                                                                                          %
%       Added conversion to [CA] for AIF with quadratic relationship                                                       %
% - LCB 7/1/2016                                                                                                           %
%       Added determineTimePoints function for brainmask to ensure consistency                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function AIF = AutoAIF_Brain(vol, paramsTime, flip, DCEflag, ConcFlag,brainmask,ss_tp,gd_tp,pk_tp)
%% Read in volume and determine matrix size
if length(size(vol)) < 4 || length(size(vol)) > 6
    warning('Incorrect matrix size!');
    return;
elseif length(size(vol)) == 5
    [nx,ny,nz,ne,nt] = size(vol);
    warning('This code only uses information for 2 echoes.'); %you'll need to code if you want to add other echoes
    display(sprintf('\nThis is a multi echo dataset with %d echoes: ', ne));
    display(sprintf('nx: %d ny: %d nz: %d nt: %d \n', nx, ny, nz, nt));
    echoFlag = 1; %flag multiecho
    TE1 = paramsTime(2); TE2 = paramsTime(3); TR = paramsTime(1);
    I = squeeze(vol(:,:,:,1,:)); %only need first echo dataset for this section
    I_v = reshape(I, [nx*ny*nz nt]);
elseif length(size(vol)) == 4
    [nx,ny,nz,nt] = size(vol);
    display(sprintf('\nThis is a single echo dataset: '));
    display(sprintf('nx: %d ny: %d nz: %d nt: %d \n', nx, ny, nz, nt));
    echoFlag = 0; %flag single echo
    TE1 = paramsTime(2); TR = paramsTime(1);
    I = vol;
    I_v = reshape(I, [nx*ny*nz nt]);
end

%% First mask brain and determine steady state and contrast arrival timepoints

%LCB added determineTimePoints function 7/1/2016
%NS deleted determineTimePoints, WHY do it twice??
ss_tp=ss_tp;
gd_tp=gd_tp;
pk_tp=pk_tp;
brainmask=brainmask;
er_brainmask = zeros(nx,ny,nz);
for z = 1:nz
    er_brainmask(:,:,z) = imerode(imfill(brainmask(:,:,z)),strel('disk', 2, 0)); %AMS changed strel('disk',5,0) to strel('disk',2,0), bigger mask
end
%figure, montage(permute(er_brainmask, [1 2 4 3])); title('Brain mask');
brainmask = er_brainmask;

peak_frame = I(:,:,:,pk_tp);

%% Calculate delta R2 star and delta R1
switch echoFlag
    case 1
        I = abs(vol(:,:,:,:,:));
        I_v = squeeze(reshape(I, [nx*ny*nz ne nt]));
        
        %calculate delta R2 star curves
        S0_TE1 = squeeze(mean(I_v(:,1,ss_tp:gd_tp),3)); 
        S0_TE1 = repmat(S0_TE1, [1 nt]);
        
        S0_TE2 = squeeze(mean(I_v(:,2,ss_tp:gd_tp),3)); 
        S0_TE2 = repmat(S0_TE2, [1 nt]);
        
        dR2s_TE1 = log(S0_TE1./squeeze(I_v(:,1,:)))/TE1; 
        dR2s_TE1(isinf(dR2s_TE1)) = 0; 
        dR2s_TE1(isnan(dR2s_TE1)) = 0;
        dR2s_TE2 = log(S0_TE2./squeeze(I_v(:,2,:)))/TE2; 
        dR2s_TE2(isinf(dR2s_TE2)) = 0; 
        dR2s_TE2(isnan(dR2s_TE2)) = 0;
        dR2s_all = (1/(TE2-TE1)).*(log( (squeeze(I_v(:,1,:))./S0_TE1) ./ (squeeze(I_v(:,2,:))./S0_TE2) ));
        dR2s_all(isinf(dR2s_all)) = 0; 
        dR2s_all(isnan(dR2s_all)) = 0;
        
        dR2s_TE1 = dR2s_TE1 .* repmat(brainmask(:), [1 nt]);
        dR2s_TE2 = dR2s_TE2 .* repmat(brainmask(:), [1 nt]);
        dR2s_all = dR2s_all .* repmat(brainmask(:), [1 nt]);
        
        clear S0_TE1 S0_TE2
        
        switch DCEflag
            case 1 %Schmainda et al. method
                % raw signal time courses for each echo
                STE1_t = squeeze(vol(:,:,:,1,:)); 
                STE2_t = squeeze(vol(:,:,:,2,:)); 
                STE1_t = reshape(STE1_t, [nx*ny*nz nt]); %vectorize
                STE2_t = reshape(STE2_t, [nx*ny*nz nt]); %vectorize
                STE1_0 = STE1_t(:,1); % t = 0 here to estimate S0 below
                STE2_0 = STE2_t(:,1); % t = 0 here to estimate S0 below
          
                % estimate of the equalibrium magnetization
                S0 = (STE1_0/sind(flip)) .* exp((TE1/(TE2-TE1)) .* log(STE1_0./STE2_0));
                S0(S0 < 1000) = 0;
                S0_t = repmat(S0, [1 nt]); %replicated for all time points
                
                % estimate 1/T2* to correct for T1 leakage
                T2star_inv = (1/(TE2-TE1)) .* log(STE1_t./STE2_t);
                
                % corrected first echo signal STE1c_t by extrapolating back to TE=0 
                STE1c_t = STE1_t .* exp(TE1*T2star_inv); 
                SBc = repmat(mean(STE1c_t(:,ss_tp:gd_tp),2), [1 nt]); %average baseline signal of corrected signal

                % DCE - delta R1 star corrected for T2/T2* and T1 effects
                clear a b
                a = S0_t*sind(flip) - STE1c_t;
                b = S0_t*sind(flip) - STE1c_t*cosd(flip);
                c = S0_t*sind(flip) - SBc*cosd(flip);
                d = S0_t*sind(flip) - SBc;
                
                deltaR1_t = (-1/TR) .* log((a./b).*(c./d));
                deltaR1_t(isinf(deltaR1_t)) = 0;
                deltaR1_t = abs(deltaR1_t);
                
                deltaR1_t = deltaR1_t .* repmat(brainmask(:), [1 nt]);

                clear a b c d
                
            case 2                 
                warning('You need to code this if you want to use T1 maps option.');
                %need to write code here for R1 curves based on T1 mapping
        end
        
    otherwise %only single echo acquisition
        I_v = abs(I_v);
        S0_TE1 = squeeze(mean(I_v(:,ss_tp:gd_tp),2)); S0_TE1 = repmat(S0_TE1, [1 nt]);
        dR2s_TE1 = log(S0_TE1./squeeze(I_v(:,:)))/TE1;
        dR2s_TE1(isinf(dR2s_TE1)) = 0; dR2s_TE1(isnan(dR2s_TE1)) = 0;
        dR2s_TE1 = dR2s_TE1 .* repmat(brainmask(:), [1 nt]);
        
        clear S0_TE1
end

%% Calculate AIF mask based on delta R2 star curves
switch echoFlag
    case 1 %multiecho
        
        %Easy thresholding
        [VX1, TX1] = max(permute(squeeze(dR2s_TE1(:,ss_tp:end)), [2 1]));
        [VX2, TX2] = max(permute(squeeze(dR2s_TE2(:,ss_tp:end)), [2 1]));
        
        mask = brainmask(:);
        mask(abs(TX1 - TX2) > 0) = 0; %get rid of noisey curves if time of peak dR2s is not the same for both echoes
        mask(VX1 > VX2) = 0; %get rid of T1 leakage effects if peak dR2s for TE2 is greater than TE1
        
        finalmask = zeros(size(mask));
        for x = 1:length(mask)
            if mask(x,1) > 0 && corr(dR2s_TE1(x,ss_tp:end)', dR2s_TE2(x,ss_tp:end)') > 0.9
                finalmask(x,1) = 1;
            end
        end
        
        loc = find(finalmask == 1);
        VXall = max(permute(squeeze(dR2s_all(loc,ss_tp:end)), [2 1]));
        loc(VXall > 100) = 0; %get rid of outliers with extreme height
        loc(loc == 0) = [];
        
    case 0 %single echo
        
        %Quick and easy time of arrival threshold
% %         [~, TX] = min(permute(squeeze(I_v(:,ss_tp:end)), [2 1])); %comment out AMS 20170911
% %         TX(TX < ss_tp) = 0; %this doesn't make sense.... %comment out AMS 20170911
% %         TX(TX > gd_tp) = 0; %comment out AMS 20170911
        [~, TX] = min(permute(squeeze(I_v(:,:)), [2 1]));
        TX(TX < gd_tp) = 0; %want min signal to occur after prebolus phase
        TX(TX >= pk_tp) = 0; %but before peak over whole brain
        TX = TX' .* brainmask(:);
        loc = find(TX > 0);
        VXTE1 = max(permute(squeeze(dR2s_TE1(loc,ss_tp:end)), [2 1]));
        loc(VXTE1 > 100) = 0; %get rid of outliers with extreme height
        loc(loc == 0) = [];
              
end

switch echoFlag
    case 1
        I_v = squeeze(I(:,:,:,2,:));
        I_v = reshape(I_v, [nx*ny*nz nt]);
end

sum_dR2 = zeros(nx*ny*nz,1);
precontrast_mean = mean(I_v(:,ss_tp:gd_tp),2);
precontrast_std = std(I_v(:,ss_tp:gd_tp), [], 2);
for xx = 1:length(loc);
    ind = loc(xx);
% %     figure(100);plot(I_v(loc(xx),:));pause(0.05);
    bolus_pts = find(squeeze(I_v(ind,1:gd_tp+10)) < precontrast_mean(ind,:)-10*precontrast_std(ind,:));
    %bolus_pts = bolus_pts + (ss_tp-1);
    if ~isempty(bolus_pts) && length(bolus_pts)*TR < 8 %AIF shouldn't enhance longer than 8 seconds - removes enhancing voxels
        base_pts = find(squeeze(I_v(ind,1:bolus_pts(1))) > precontrast_mean(ind,:)-3*precontrast_std(ind,:));
        if ~isempty(base_pts)
            arrival_tmp = bolus_pts(1) - base_pts(end);
            if arrival_tmp <= 3 && arrival_tmp > 0
                switch echoFlag
                    case 1
                        sum_dR2(ind,:) = trapz(dR2s_all(ind,base_pts(end):bolus_pts(end)));
                    case 0
                        sum_dR2(ind,:) = trapz(dR2s_TE1(ind,base_pts(end):bolus_pts(end)));
                end
            end
            clear arrival_tmp
        end
    end
end

new_loc = find(sum_dR2 > 0);
table = cat(2, sum_dR2(new_loc,:), new_loc);
table = abs(sortrows(-table, 1));
switch echoFlag
    case 1
        final_loc = table(1:min(10,size(table,1)),2); %take top 10 -- EDIT AMS20160216 - if <10 are found, take min of (10,#found)
        if size(table,1)<10 % ADDED AMS20160216 - notify user of # of voxels found
            display(sprintf('\nFound %d potential AIF voxels, keeping all voxels.', size(table,1)));
            warning('Found less than 10 potential AIF voxels, proceeding with AIF calculation');
        else
            display(sprintf('\nFound %d potential AIF voxels, keeping top 10.', size(table,1)));
        end
    case 0
        final_loc = table(1:min(20,size(table,1)),2); %take top 10 -- EDIT AMS20160216 - if <10 are found, take min of (10,#found)
        if size(table,1)<20 % ADDED AMS20160216 - notify user of # of voxels found
            display(sprintf('\nFound %d potential AIF voxels, keeping all voxels.', size(table,1)));
            warning('Found less than 20 potential AIF voxels, proceeding with AIF calculation');
        else
            display(sprintf('\nFound %d potential AIF voxels, keeping top 20.', size(table,1)));
        end
end

    

AIFmask = zeros(nx*ny*nz,1);
AIFmask(final_loc,:) = 1;
AIFmask = reshape(AIFmask, [nx ny nz]);
save AIFmask.mat AIFmask

%% Plot AIF Curve
% switch echoFlag
%     case 1
%         h = figure; hold on;
%         title('AIF delta R2 star');
%         plot(mean(dR2s_all(final_loc,ss_tp:end)), 'k', 'LineWidth', 2);
%         plot(mean(dR2s_TE1(final_loc,ss_tp:end)), '--k', 'LineWidth', 2);
%         plot(mean(dR2s_TE2(final_loc,ss_tp:end)), ':k', 'LineWidth', 2);
%         legend('Dual Echo', 'TE1', 'TE2');
%         xlabel('Time points', 'FontSize', 14); ylabel('delta R2 star 1/second', 'FontSize', 14);     
%     case 0
%         h = figure; hold on;
%         title('AIF delta R2 star');
%         plot(mean(dR2s_TE1(final_loc,ss_tp:end)), 'k', 'LineWidth', 2);
%         xlabel('Time points', 'FontSize', 14); ylabel('delta R2 star 1/second', 'FontSize', 14);
% end
% 
% savefig(h, 'AIF_dR2s_timecourse.fig');

%For visualization
% figure('Visible', 'off'), hm = montage(permute(AIFmask, [1 2 4 3])); aif = hm.CData;  
% figure('Visible', 'off'), hm = montage(permute(peak_frame, [1 2 4 3])); peak = hm.CData; close;
% h2 = figure; imagesc(imfuse(aif, peak)); axis off;
% savefig(h2, 'AIFlocations.fig');

%% Calculate AIF mask based on delta R1 star curves

% Ashley mentioned she uses a different curve. My current dataset is not
% ideal for this...maybe add this feature in the future.

% output delta R1 curves based on the AIF mask determined by delta R2star.
% Most likely garbage based on dual echo spiral datasets. LCB 12/2015

%% Convert to concentration 
% AMS 6/20/2016 - added section to convert to [CA]

switch ConcFlag
    case 1
        switch echoFlag
            case 1
                display('Calculating AIF in concentration');
                %Compute AIF concentration curve using quadratic relationship
                %Calamante et al MRM 2009; 61:486

                CAIF=zeros(size(final_loc,1),size(dR2s_all,2));
                % % a1=7.62; a2=0.574; htcvar=1.1378*(0.4)/(1-0.4)^2; %1.5T
                a1=0.493; 
                a2=2.62; 
                htcvar=1.1378*(0.4)/(1-0.4)^2; % 3T
                for j=1:size(final_loc,1)
                     for n=1:size(dR2s_TE1,2)                
                        F=[a2*htcvar,a1,-dR2s_all(final_loc(j),n)];
                        C=roots(F);
                        if isreal(C(1)) == 1 %redundant, this may make mistakes
                         CAIF(j,n)=nanmax(C);
                        end
                     end
                end

                trend=nanmean(CAIF(:,ss_tp:gd_tp),2);
                Caifq_avg_dt=CAIF-repmat(trend,[1 nt]); 

                Caifq_DCE=Caifq_avg_dt; 

            case 0
                display('Calculating AIF in concentration');
                %Compute AIF concentration curve using quadratic relationship
                %Calamante et al MRM 2009; 61:486

                CAIF=zeros(size(final_loc,1),size(dR2s_TE1,2));
                % % a1=7.62; a2=0.574; htcvar=1.1378*(0.4)/(1-0.4)^2; %1.5T
                a1=0.493; a2=2.62; htcvar=1.1378*(0.4)/(1-0.4)^2; % 3T
                for j=1:size(final_loc,1)
                     for n=1:size(dR2s_TE1,2)                
                        F=[a2*htcvar,a1,-dR2s_TE1(final_loc(j),n)];
                        C=roots(F);
                        if isreal(C(1)) == 1
                         CAIF(j,n)=nanmax(C);
                        end
                     end
                end

                trend=nanmean(CAIF(:,ss_tp:gd_tp),2);
                Caifq_avg_dt=CAIF-repmat(trend,[1 nt]); 

                Caifq_DCE=Caifq_avg_dt; 
        end

end

%% Create ouput
switch echoFlag
    case 1
        AIF = mean(dR2s_all(final_loc,:));
        
        switch ConcFlag
            case 1
                display('AIF output in DR2s (first row) and concentration (second row)');
%                 AIF(2,:) = mean(Caifq_DCE);
                AIF(2,:) = nanmean(Caifq_avg_dt);
        end
        
        
        switch DCEflag
            case 1
                AIF = cat(1, mean(dR2s_all(final_loc,:)), mean(deltaR1_t(final_loc,:)));
            case 2
                warning('You need to code this if you want to use T1 maps option.');
        end
        
    otherwise
        AIF = mean(dR2s_TE1(final_loc,:));
        switch ConcFlag
            case 1
                display('AIF output in DR2s (first row) and concentration (second row)');
                AIF(2,:) = nanmean(Caifq_DCE);
        end
end

end

