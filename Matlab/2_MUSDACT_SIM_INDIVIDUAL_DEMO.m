%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                         VISUALISE, TRIM AND SEGMENT                          %
%                  MEAN ABSOLUTE DEVIATION AND CLASSIFICATION                  %
%                  OF WEARABLE ACCELEROMETRY FOR INDIVIDUALS                   %
%                                                                              %
%                               25 OCTOBER, 2022                               %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CITATION A.P.A

% Mendoza, J.I., Danso, A., Luck, G., Rantalainen, T., Palmberg., L, & Chastin, S. (2022). 
% Musification of Accelerometry Data Towards Raising Awareness of Physical Activity. 
% Conference on Sonification of Health and Environmental Data SoniHED, 
% 27-28 October 2022, Stockholm, Sweden.

% ------------------------------------------------------------------------------

% INFORMATION

% Dependencies:
%     novelty
%     peaksfind

% Tested with Matlab R2017b and R2019b.

% Instructions:
%     Edit the values indicated with an arrow like this: <---
%     Run the program, close your eyes and hope for the best.

% ==============================================================================
% Initialisation:

clc
clear
close all
restoredefaultpath

% ------------------------------------------------------------------------------
% Declare paths:

   info.paths.mcode = ''; % <--- folder where this program and author's functions are
   info.paths.data  = ''; % <--- path for data (raw and processed)

addpath(genpath(info.paths.mcode))
cd(info.paths.data)
% =============================================================================
clc 
close all
i_p = 0;

% ------------------------------------------------------------------------------
% Input data information and parameters:

% ..............................................................................

common.day_start  = 0; % <---
common.hour_start = 0; % <---
common.day_end    = 1; % <---
common.hour_end   = 0; % <---

% ..............................................................................
% Participants' data:

i_p = i_p + 1;
info.data(i_p).fname_in    = 'SIM_LOW.mat';   % <--- file with MAD and activity classes, with extension
info.data(i_p).descr.short = 'LO';           % <---
info.data(i_p).descr.long  = 'LOW ACTIVITY'; % <---
info.data(i_p).mn_offset   = 0 ; % <--- offset in 5s epochs (h.m.s) so that day 1 and hour 0 is 00:00:00 (H:MM:SS)

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

i_p = i_p + 1;
info.data(i_p).fname_in    = 'SIM_HIGH.mat';    % <--- file with MAD and activity classes, with extension
info.data(i_p).descr.short = 'HI';            % <---
info.data(i_p).descr.long  = 'HIGH ACTIVITY'; % <---
info.data(i_p).mn_offset   = 0 ; % <--- offset in 5s epochs (h.m.s) so that day 1 and hour 0 is 00:00:00 (H:MM:SS)

% ..............................................................................


info.activity_lbl = {'lying','sitting','upright','walking'};  % <--- class labels (0:3)

% ..............................................................................

process_participants = [1,2]; % <--- process only these participants

info.int_param = []; % <--- integrate with sliding window and compute log: [ width , hop ] in 5s. epochs; [] = don't

% info.int_param = [ 12 , 6 ];  % [ 1 , 0.5 ] minute 
% info.int_param = [ 60 , 30 ]; % [ 5 , 2.5 ] minutes
info.int_param = [ 120 , 60 ]; %  [ 10 , 5  ] minutes [[[[[[[[[[ WORKS WELL WITH 24-HOUR, 3-VOICE SONIFICATION, 30-90 s. ]]]]]]]]]]

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

% The following parameters apply only if data is integrated (info.int_param not empty):

% info.seg.bw          = 2:4:64; % <--- novelty score bandwidth (*) in time units (**)
info.seg.bw          = 2:2:32; % <--- novelty score bandwidth (*) in time units (**) [[[[[[[[[[ WORKS WELL WITH 24-HOUR, 3-VOICE SONIFICATION, 30-90 s. ]]]]]]]]]]
% info.seg.bw          = 2:2:16; % <--- novelty score bandwidth (*) in time units (**)
info.seg.theta       = 0; % <--- novelty threshold for segmentation boundaries
info.seg.delta_d     = 2; % <--- tolerance to rectify drifting segmentation boundaries in time units (**)
info.seg.f_corr_marg = 3; % <--- fraction of the coarser novelty bandwidth to search from the edges for rogue segmentation boundaries
info.seg.max_gran    = 3; % <--- maximum number of granularity levels allowed (finest and coarser are always retained)
info.seg.seg_func    = 'median'; % <---- segments' Activity Score ('sum','mean','median')

% (*) bandwidth is the standard deviation of the gaussian that tapers the checkerboard kernel
% (**) time units is window lengths if integrated, otherwise epochs of 5s.

vis.make.a_score = 1; % <--- display activity score (***); 0 = don't
vis.make.classes = 1; % <--- display integrated classes; 0 = don't
vis.make.simmat  = 1; % <--- display similarity matrix; 0 = don't
vis.make.nov     = 1; % <--- display novelty scores; 0 = don't
vis.make.pk_raw  = 1; % <--- display raw novelty peaks (drifting segmentation boundaries); 0 = don't
vis.make.sb_rec  = 1; % <--- display rectified segmentation boundaries; 0 = don't
vis.make.seg_act = 1; % <--- display segmented activity; 0 = don't

vis.opt.i_mon    = 2;  % <--- monitor to display figures
vis.opt.slow_rec = 0; % <--- 1 = slowly display rectification of segmentation boundaries and elimination of rogue boundaries; 0 = don't
vis.opt.gray     = 0; % <--- 1 = display in grayscale; 0 = use default chromatic colours

% (***) Activity Score = log3( 1 + int(MAD) )

write_gif_size = [];
% write_gif_size = [76,1240]; % <--- [height,width] = position of GIF file to save of segmented activity heatmap; [] = don't save GIF
     write_csv = 1;          % <--- save space-separated CSV files

% ------------------------------------------------------------------------------
if info.seg.max_gran < 2
    warning('info.seg.max_gran should be 2 or more.')
    return
end

 info.dh_lbl = sprintf( 'd[%i_%i]_h[%i_%i]' , ...
                         common.day_start   , ...
                         common.day_end     , ...
                         common.hour_start  , ...
                         common.hour_end      ...
                         );

n_class = length(info.activity_lbl);
n_proc_p = length(process_participants);
n_bw = length(info.seg.bw);

out = struct;

monpos = get(0, 'MonitorPositions');
scrsz = monpos(vis.opt.i_mon,:);

fig_w = scrsz(3)/n_proc_p;
fig_pos = scrsz;
fig_pos(3) = fig_w;

if isempty(info.int_param)
    
    n_sp = 2;
else
    
    fldn = fieldnames(vis.make);
    n_sp = 0;
    for i_fldn = 1:length(fldn)
        
        if vis.make.(fldn{i_fldn}) 
            if strcmp( fldn{i_fldn} , 'classes')
                n_sp = n_sp + n_class;
            else
                n_sp = n_sp + 1;
            end
        end
    end
    
    if n_sp
        
        h_sp{n_sp} = [];
    end
end

if vis.opt.gray && (n_sp || isempty(info.int_param))
    
    gray_cmap.simmat = repmat( flip(0.5:0.00196:1)' ,1,3);
    gray_cmap.nov = repmat( (flip(0.9:0.000392:1)').^20,1,3);
    gray_cmap.sb = repmat( [1;0.8;0],1,3);
    gray_cmap.seg = repmat( rescale( asin(flip(0:0.001565:1)'), 0, 1) ,1,3);
end

fig_1_h{n_proc_p} = [];

for i_p = process_participants

    if ~isfield(info.data(i_p),'loaded') || isempty(info.data(2).loaded)
      
        load(info.data(i_p).fname_in);
        info.data(i_p).activities = activities;
        info.data(i_p).hours_dec = double( (activities.tStamps-activities.tStamps(1))/1000 ) / (60*60); % activities.tStamps are 5-second epochs
        info.data(i_p).loaded = 1;
        disp('data loaded')
    end
    [~,this_fname] = fileparts(info.data(i_p).fname_in);

     day_start = common.day_start;
    hour_start = common.hour_start;
       day_end = common.day_end;
      hour_end = common.hour_end;
      
     mn_offset = info.data(i_p).mn_offset;

    i_start_seg = mn_offset + day_start*17280 + hour_start*720 + 1;
    i_end_seg = mn_offset + day_end*17280 + hour_end*(720);
    i_selection = (i_start_seg:i_end_seg);
    
    out(i_p).noint.sel_mad = info.data(i_p).activities.ukk.mad(i_selection); % time-selected mean absolute deviation
    out(i_p).noint.sel_class = info.data(i_p).activities.classification(i_selection) + 1; % time-selected class, classes start at 0, they are shifted so first class = 1
    out(i_p).noint.sel_class( isnan(out(i_p).noint.sel_class) ) = 0; % turn NaN into zeros
    
    start_datetime = datetime( info.data(i_p).activities.tStamps(i_start_seg) / 1000 ,'convertfrom','posixtime','Format','HH:mm:ss');
    start_datetime_vec = datevec( start_datetime );
    start_time_vec = start_datetime_vec(4:end);
    end_datetime = datetime( info.data(i_p).activities.tStamps(i_end_seg) / 1000 ,'convertfrom','posixtime','Format','HH:mm:ss');
    fprintf('start/end times: %s / %s \n',char(start_datetime),char(end_datetime))
    
    length_sel = length(i_selection);
    i_adjusted = (1:length_sel) + hour_start*720 ;
    out(i_p).noint.sel_hours_dec = info.data(i_p).hours_dec( i_adjusted ); % 5-second epochs
    
    if isempty(info.int_param)

        plot_sel_hours_dec = out(i_p).noint.sel_hours_dec;
        plot_sel_mad = out(i_p).noint.sel_mad;
    else 
        
        % ......................................................................
        % Integrate silding window:
        
        n_windows = fix((length_sel - info.int_param(1)) / info.int_param(2) + 1);
        half_w = info.int_param(1) / 2;

        out(i_p).int.sel_hours_dec = zeros(n_windows,1);
        out(i_p).int.sel_mad = zeros(n_windows,1);
        out(i_p).int.sel_class = zeros(n_windows,n_class);
        
        bool_sel_class = zeros(length_sel,n_class);
        
        for i_class = 1:n_class
           
            bool_sel_class(:,i_class) = out(i_p).noint.sel_class == i_class;
        end

        for i_w = 1:n_windows

            i_this_window = ( info.int_param(2) * (i_w - 1) ) + ( 1 : + info.int_param(1) );
                        
            out(i_p).int.sel_hours_dec(i_w) = out(i_p).noint.sel_hours_dec( i_this_window(end) - half_w ); % new decimal hours at half of the window

            out(i_p).int.sel_mad(i_w) = sum( out(i_p).noint.sel_mad( i_this_window ) ); % integrate MAD
            
            for i_class = 1:n_class
                
                out(i_p).int.sel_class(i_w,i_class) = sum( bool_sel_class( i_this_window , i_class ) ); % integrate class (activity)
            end
        end

        out(i_p).int.act_score = log( 1 + out(i_p).int.sel_mad) / log(3); % Activity Score
                
        plot_sel_hours_dec = out(i_p).int.sel_hours_dec;
        plot_sel_mad = out(i_p).int.act_score;

        % ......................................................................
    end
    
    i_sp = 1;
    
    if n_sp || isempty(info.int_param)
        
        fig_pos(1) = fig_pos(1) + (fig_w*(i_p-1));
        fig_1_h{i_p} = figure('Position',fig_pos);
        pause(0.8) % increase this value if selected monitor is 2 but displays in monitor 1
        fig_1_h{i_p}.Position = fig_pos;
                        
        subplot(n_sp,1,i_sp)
        if vis.opt.gray
            h_sp{i_sp} = plot(plot_sel_hours_dec,plot_sel_mad,'k','LineWidth',2);
        else
            h_sp{i_sp} = plot(plot_sel_hours_dec,plot_sel_mad,'LineWidth',2);
        end
        xlim([plot_sel_hours_dec(1),plot_sel_hours_dec(end)])
        hour_step = diff(h_sp{i_sp}.Parent.XTick(1:2));
        old_XTick = h_sp{i_sp}.Parent.XTick;
        new_XTick_lbl = old_XTick;
        i_day_overflow = new_XTick_lbl > 24;
        while any(i_day_overflow)
            hour_start = new_XTick_lbl( find(i_day_overflow,1) ) - 24;
            overflow_replacement = hour_start:hour_step:(hour_step * sum(i_day_overflow));
            new_XTick_lbl(i_day_overflow) = overflow_replacement;
            i_day_overflow = new_XTick_lbl > 24;
        end
        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
    end
    
    if vis.make.a_score
    
        xlabel('hour')
        if isempty(info.int_param)
            title(sprintf('MAD for ''%s'' (%s)',this_fname,info.data(i_p).descr.short));
        else
            title(sprintf('Activity Score for ''%s'' (%s)',this_fname,info.data(i_p).descr.short));
        end
    else
        i_sp = i_sp - 1;
    end

    if isempty(info.int_param)
        
        i_sp = i_sp + 1;
        subplot(n_sp,1,2)
        if vis.opt.gray
            h_sp{i_sp} = plot(plot_sel_hours_dec,out(i_p).noint.sel_class,'k.','MarkerSize',25);
        else
            h_sp{i_sp} = plot(plot_sel_hours_dec,out(i_p).noint.sel_class,'.','MarkerSize',25);
        end
        xlim([plot_sel_hours_dec(1),plot_sel_hours_dec(end)])
        ylim([1, max(out(i_p).noint.sel_class)]) % exclude zeros (formerly NaN) from axis
        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
        these_class = unique(out(i_p).noint.sel_class);
        these_class( these_class == 0) = [];
        h_sp{i_sp}.Parent.YTick = these_class;
        h_sp{i_sp}.Parent.YTickLabel = info.activity_lbl(these_class);
        xlabel('hour')
        title('Activity States');
    else

        if vis.make.classes
    
            for i_class = 1:n_class
                
                i_sp = i_sp + 1;
                subplot(n_sp,1,i_sp)
                if vis.opt.gray
                    h_sp{i_sp} = plot(plot_sel_hours_dec,out(i_p).int.sel_class(:,i_class),'k','LineWidth',2); 
                else
                    h_sp{i_sp} = plot(plot_sel_hours_dec,out(i_p).int.sel_class(:,i_class),'LineWidth',2); 
                end
                xlim([out(i_p).int.sel_hours_dec(1),out(i_p).int.sel_hours_dec(end)])
                h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
                xlabel('hour')
                title(info.activity_lbl{i_class});
            end
        end
    end

    if ~isempty(info.int_param) 
        
        % ......................................................................
        % Multigranular temporal segmentation boundaries:
        
        out(i_p).int.nov = zeros(n_bw,n_windows);
        out(i_p).int.sim_mat = squareform(pdist(out(i_p).int.sel_class,'cityblock'));
        out(i_p).int.i_seg.all =  zeros(n_bw,n_windows);
        
        for i_bw = 1:n_bw
            
            out(i_p).int.nov(i_bw,:) = novelty(out(i_p).int.sim_mat,info.seg.bw(i_bw),0);
            out(i_p).int.nov(i_bw,:) = conv( out(i_p).int.nov(i_bw,:) , gausswin(info.seg.bw(i_bw)) ,'same');
            i_pk = peaksfind( out(i_p).int.nov(i_bw,:) ,info.seg.theta);
            out(i_p).int.i_seg.each{i_bw} = i_pk;
            out(i_p).int.i_seg.all(i_bw,i_pk) = 1;
        end
        
        if n_sp
            
            i_map_xticks = zeros(1,length(old_XTick));
            xticks_dist = (abs(plot_sel_hours_dec - old_XTick));
            
            for i_mxt = 1:length(old_XTick)
                
                i_map_xticks(i_mxt) = find( xticks_dist(:,i_mxt) == min(xticks_dist(:,i_mxt)) );
            end
        end

        if vis.make.simmat    
            i_sp = i_sp + 1;
            subplot(n_sp,1,i_sp)
            h_sp{i_sp} = imagesc(out(i_p).int.sim_mat);
            h_sp{i_sp}.Parent.YTick = [];
            h_sp{i_sp}.Parent.XTick = i_map_xticks;
            h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
            xlabel('hour')
            title('Activity Self-Similarity')
            if vis.opt.gray
                colormap(h_sp{i_sp}.Parent,gray_cmap.simmat)
            end
        end
        
        if vis.make.nov
            i_sp = i_sp + 1;
            subplot(n_sp,1,i_sp)
            h_sp{i_sp} = imagesc(out(i_p).int.nov);
            h_sp{i_sp}.Parent.XTick = i_map_xticks;
            h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
            ylabel('granularity')
            xlabel('hour')
            title('Activity Novelty')
            if vis.opt.gray
                colormap(h_sp{i_sp}.Parent,gray_cmap.nov)
            end
        end
        
        if vis.make.pk_raw
            i_sp = i_sp + 1;
            subplot(n_sp,1,i_sp)
            h_sp{i_sp} = imagesc(out(i_p).int.i_seg.all);
            h_sp{i_sp}.Parent.XTick = i_map_xticks;
            h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
            ylabel('granularity')
            xlabel('hour')
            title('Activity Boundaries')
            if vis.opt.gray
                colormap(h_sp{i_sp}.Parent,gray_cmap.sb)
            end
        end
        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        % Rectify drifting boundaries:

        out(i_p).int.i_seg.each_rect = out(i_p).int.i_seg.each; % vectors "bw" of indexes "pk" in matrix "all"
        out(i_p).int.i_seg.all_rect = out(i_p).int.i_seg.all; % matrix where raw boundaries = 1
        out(i_p).int.i_seg.all_rect(1,out(i_p).int.i_seg.each_rect{1}) = 2; % update matrix: mark finest granularity boundaries as rectified (=2)
        
        if vis.make.sb_rec
            
            i_sp = i_sp + 1;
        end

        for i_sb_1 = 1:length(out(i_p).int.i_seg.each_rect{1}) % values in first vector
            
            i_pk_1 = out(i_p).int.i_seg.each_rect{1}(i_sb_1);
            i_sb_n = i_sb_1;
            i_bw = 1;
            
            while i_bw < n_bw
                
                i_pk_prev = out(i_p).int.i_seg.each{i_bw}(i_sb_n);
                
                if isempty(i_pk_prev)
                    break
                end
                
                i_bw = i_bw + 1;
                
                i_sb_n = (  (out(i_p).int.i_seg.each{i_bw} >= (i_pk_prev-info.seg.delta_d) ) & (out(i_p).int.i_seg.each{i_bw} <= (i_pk_prev+info.seg.delta_d) ));
                
                if any(i_sb_n)
                    
                    smaller_info.seg.delta_d = info.seg.delta_d;
                    
                    while sum(i_sb_n) > 1 % if boundaries are closer than tolerance margin
                        
                        smaller_info.seg.delta_d = smaller_info.seg.delta_d - 1;
                        i_sb_n = (  (out(i_p).int.i_seg.each{i_bw} >= (i_pk_prev-smaller_info.seg.delta_d) ) & (out(i_p).int.i_seg.each{i_bw} <= (i_pk_prev+smaller_info.seg.delta_d) ));
                    end
                    
                    i_pk_n = out(i_p).int.i_seg.each_rect{i_bw}(i_sb_n); % non-rectified value in vector
                    out(i_p).int.i_seg.all_rect(i_bw,i_pk_n) = 0; % update matrix: remove non-rectified boundary
                    
                    out(i_p).int.i_seg.each_rect{i_bw}(i_sb_n) = i_pk_1; % replace (rectify) value in vector with finest granularity value
                    out(i_p).int.i_seg.all_rect(i_bw,i_pk_1) = 2; % update matrix: insert rectified boundary
                    
                    if vis.make.sb_rec
                        
                        subplot(n_sp,1,i_sp)
                        h_sp{i_sp} = imagesc( out(i_p).int.i_seg.all_rect );
                        h_sp{i_sp}.Parent.XTick = i_map_xticks;
                        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
                        if vis.opt.slow_rec
                            if vis.opt.gray
                                colormap(h_sp{i_sp}.Parent,gray_cmap.sb)
                            end
                            drawnow
                        end
                    end
                else
                    break
                end
            end
        end
        
        out(i_p).int.i_seg.all_rect( out(i_p).int.i_seg.all_rect == 1) = 0; % remove unrectified boundaries (false boundaries at the edges)
        out(i_p).int.i_seg.all_rect( out(i_p).int.i_seg.all_rect == 2) = 1; % mark all rectified boundaries as 1
        
        if vis.make.sb_rec
            
            if vis.opt.slow_rec
                pause(1)
                drawnow
            end
            subplot(n_sp,1,i_sp)
            h_sp{i_sp} = imagesc( out(i_p).int.i_seg.all_rect );
            h_sp{i_sp}.Parent.XTick = i_map_xticks;
            h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
            if vis.opt.gray
                colormap(h_sp{i_sp}.Parent,gray_cmap.sb)
            end
        end
        
        for i_bw = 2:n_bw
            
            out(i_p).int.i_seg.each_rect{i_bw} =  find(out(i_p).int.i_seg.all_rect(i_bw,:)); % update vectors
        end
        
        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        % Eliminate redundancy and reduce to maximum number of granularities:
        
        max_n_gran = info.seg.max_gran - 2;
        n_b = sum( out(i_p).int.i_seg.all_rect,2 );
        diff_n_b =  abs(diff( n_b));
        [d_sorted,i_d_sort] = sort(diff_n_b);
        i_start_nonzero = find(d_sorted,1); % zero diffs are redundant
        if max_n_gran >= n_bw
            max_n_gran = n_bw - 1;
        end
        i_redux = unique([ 1; sort( i_d_sort(i_start_nonzero : i_start_nonzero + max_n_gran - 1) ) ; n_bw ]);
        out(i_p).int.i_seg.all_rect_rdx = out(i_p).int.i_seg.all_rect( i_redux , : );
        
        if vis.make.sb_rec
            
            if vis.opt.slow_rec
                pause(0.5)
                drawnow
            end
            subplot(n_sp,1,i_sp)
            h_sp{i_sp} = imagesc( out(i_p).int.i_seg.all_rect_rdx );
            h_sp{i_sp}.Parent.XTick = i_map_xticks;
            h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
            ylabel('granularity')
            xlabel('hour')
            title('Rectified and Reduced Activity Boundaries')
            if vis.opt.gray
                colormap(h_sp{i_sp}.Parent,gray_cmap.sb)
            end
        end

        out(i_p).int.i_seg.each_rect_rdx{1} = out(i_p).int.i_seg.each_rect{1};
        for i_bw = 2:length(i_redux)
            
            out(i_p).int.i_seg.each_rect_rdx{i_bw} = find(out(i_p).int.i_seg.all_rect_rdx(i_bw,:)); % make vector cells of rectified boundaries
        end
        
        % ......................................................................
        % Segmented activity:
        
        n_gran = length(out(i_p).int.i_seg.each_rect_rdx);
        out(i_p).int.seg_act.all = zeros( size( out(i_p).int.i_seg.all_rect_rdx )); % initialise heatmap

        for i_gran = 1:n_gran
            
            bounds_n_ending = [ out(i_p).int.i_seg.each_rect_rdx{i_gran} , n_windows ];
            i_pk_end = 0;
            
            for i_sb = 1:length( bounds_n_ending )
                
                i_pk_start = i_pk_end + 1;
                i_pk_end = bounds_n_ending(i_sb);
                
                if strcmp(info.seg.seg_func,'sum')
                    
                    out(i_p).int.seg_act.each{i_gran}(i_sb) = sum( out(i_p).int.act_score(i_pk_start:i_pk_end) );
                    
                elseif strcmp(info.seg.seg_func,'mean')
                    
                    out(i_p).int.seg_act.each{i_gran}(i_sb) = mean( out(i_p).int.act_score(i_pk_start:i_pk_end) );
                    
                elseif strcmp(info.seg.seg_func,'median')
                    
                    out(i_p).int.seg_act.each{i_gran}(i_sb) = median( out(i_p).int.act_score(i_pk_start:i_pk_end) );
                end
                
                out(i_p).int.seg_act.all(i_gran,i_pk_start:i_pk_end) = out(i_p).int.seg_act.each{i_gran}(i_sb); % update matrix
            end
        end
        
        if vis.make.seg_act
            
            i_sp = i_sp + 1;
            subplot(n_sp,1,i_sp)
            h_sp{i_sp} = imagesc( out(i_p).int.seg_act.all );
            h_sp{i_sp}.Parent.XTick = i_map_xticks;
            h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
            ylabel('granularity')
            xlabel('hour')
            title('Segmented Activity Score')
            if vis.opt.gray
                colormap(h_sp{i_sp}.Parent,gray_cmap.seg)                
            end
        end 
        
        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        % Make sequencing matrix:
        
        % cols = [ st_1 , d_1, a_1, d_2, a_2, d_g, a_g , ... , d_n, a_n ]
        
        % st_1 = time (windows) of the start of the segment at finest granularity (g = 1)
        %  d_g = duration (windows) of the segment at granularity g
        %  a_g = integrated activity of the segment at granularity g
        
        % If a new segment doesn't start, then d_g and a_g are NaN
        
        out(i_p).int.sequence = NaN( length( out(i_p).int.i_seg.each_rect_rdx{1} ) + 1, 2 * n_gran + 1); % initialise sequence matrix
        
        out(i_p).int.sequence(:,1) = [ 0, out(i_p).int.i_seg.each_rect_rdx{1} ]' ;
        out(i_p).int.sequence(:,2) = diff( [ out(i_p).int.sequence(:,1) ; n_windows ] );
        out(i_p).int.sequence(:,3) = out(i_p).int.seg_act.each{1}';
        
        for i_st_1 = 1:( length( out(i_p).int.i_seg.each_rect_rdx{1} ) + 1)
            
            i_start_seg = out(i_p).int.sequence(i_st_1,1);
            
            i_col = 3;
            
            for i_gran = 2:n_gran
                
                % duration:
                i_col = i_col + 1;
                
                if i_st_1 == 1
                    
                    out(i_p).int.sequence(i_st_1,i_col) = diff([ 0 , out(i_p).int.i_seg.each_rect_rdx{i_gran}(1) ]);
                else
                    
                    if  out(i_p).int.i_seg.all_rect_rdx( i_gran, i_start_seg)
                        
                        bounds_this_g = [ out(i_p).int.i_seg.all_rect_rdx(i_gran,:) , n_windows ];
                        
                        i_end_seg = i_start_seg + find( bounds_this_g(i_start_seg + 1: end) , 1);
                        out(i_p).int.sequence(i_st_1,i_col) = diff( [ i_start_seg , i_end_seg ] );
                    else
                        break
                    end
                end
                
                % int. act.:
                i_col = i_col + 1;
                
                if i_st_1 == 1
                    
                    out(i_p).int.sequence(i_st_1,i_col) = out(i_p).int.seg_act.each{i_gran}(1);
                else
                    
                    i_seg_g = find( out(i_p).int.i_seg.each_rect_rdx{i_gran} == i_start_seg ) + 1;
                    out(i_p).int.sequence(i_st_1,i_col) = out(i_p).int.seg_act.each{i_gran}(i_seg_g);
                end
            end
        end
        % ......................................................................
            
        if write_csv || ~isempty(write_gif_size) % make file name addendum
            
            fname_add = sprintf( '_%s_%s'                  , ...
                                info.data(i_p).descr.short , ...
                                info.dh_lbl                  ...
                                );
        end
        
        if write_csv
            
            % space-delimited CSV single-precision floats:
            full_fn_csv = [info.paths.data,'/',this_fname,fname_add,'.csv'];
            fid = fopen(full_fn_csv,'wt');
            
            % header:
            fprintf( fid , '%i %i %.7g %i\n' , ...
                     n_windows                                 , ...
                     ceil( sum( out(i_p).int.seg_act.all(:) )) , ...
                     mean( out(i_p).int.seg_act.all(:) )       , ...
                     n_gran                                      ...
                     ); 
            
            for i_csv = 1:size(out(i_p).int.sequence,1)
                
                if isempty(info.int_param)
                    
                    fprintf(fid,'%.7g %i\n',out(i_p).noint.sel_mad(i_csv),out(i_p).noint.sel_class(i_csv));
                else
                    fprintf(fid,'%.7g ',out(i_p).int.sequence(i_csv,:));
                    fprintf(fid,'\b\n');
                end
            end
            fclose(fid);
        end

        if ~isempty(write_gif_size)
            
            fig_2_h = figure('visible','off');
            hmap = imagesc( out(i_p).int.seg_act.all );
            
            frame = getframe( hmap.Parent );
            im = frame2im(frame);
            im = imresize(im,write_gif_size);
            [imind,cm] = rgb2ind(im,256);
            
            full_fn_gif = [info.paths.data,'/',this_fname,fname_add,'.gif'];
            imwrite(imind,cm,full_fn_gif,'gif')
        end
    end
end
    
disp(' ')
disp('done')