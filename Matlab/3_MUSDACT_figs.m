%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                                 MAKE FIGURES                                 %
%                                                                              %
%                         VISUALISE, TRIM AND SEGMENT                          %
%                  MEAN ABSOLUTE DEVIATION AND CLASSIFICATION                  %
%                  OF WEARABLE ACCELEROMETRY FOR INDIVIDUALS                   %
%                                                                              %
%                                25 OCTOBER, 2022                              %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INFORMATION

% Description:
%     This program makes beautiful figures to display the data produced by    
%     "MUSDACT_SIM_INDIVIDUAL_DEMO.m"

% Dependencies:
%     MUSDACT_SIM_INDIVIDUAL_DEMO.m

% Tested with Matlab R2017b and R2019b.

% Instructions:
%     Edit the values indicated with an arrow like this: <--- (length of the arrow may vary)
%     Run the program "MUSDACT_SIM_INDIVIDUAL_DEMO.m"
%     Run this program, close your eyes and hope for the best.

% =============================================================================
clc 
close all
fig.make = zeros(9,1);
fig.opt.lbl = '';

% ------------------------------------------------------------------------------
% Input parameters:

fig.opt.savepath = ''; % <--- path to save figures

fig.opt.titles   = 0;  % <--- 1 = display full titles, 0 = only long description (*)
fig.opt.gray     = 1;  % <--- 1 = display in grayscale, 0 = use default chromatic colours
fig.opt.headroom = 10; % <--- fraction of data maxima for plotting headroom
fig.opt.ln_wdth  = 1;  % <--- line width
fig.opt.mk_size  = 6;  % <--- marker size
fig.opt.pnl_lbl  = 1;  % <--- 1 = display an alphabetic label for each subplot, 0 = don't
fig.opt.sp_size  = [361,78]; % <--- sub-plot [width,height]
fig.opt.i_mon    = 2;  % <--- monitor to display figure, monitor size is figure's size (will use 2nd monitor if bigger)
fig.opt.size_fig = [1080,1920]; % <--- [width, height] of figure with all sub-plots; [] = use monitor's size
% fig.opt.size_fig = []; 
fig.opt.dpi      = 300; % <--- resolution of TIFF file (dpi)
fig.opt.h_resc   = 0.8; % <--- horizontal rescaling factor for TIFF file

% (*) info.data(i_p).descr.long

% 1 = display, 0 = don't (comment/uncomment):
fig.make(1) = 1; % <--- MAD
fig.make(2) = 1; % <--- activities (classes)
fig.make(3) = 1; % <--- activity score (***)
fig.make(4) = 1; % <--- integrated classes
fig.make(5) = 1; % <--- similarity matrix
fig.make(6) = 1; % <--- novelty scores
fig.make(7) = 1; % <--- raw novelty peaks (drifting segmentation boundaries)
fig.make(8) = 1; % <--- rectified segmentation boundaries
fig.make(9) = 1; % <--- segmented activity

% fig.opt.lbl = 'all'; % <--- 'label' for TIFF filename, empty or commented = don't save

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  
% alternatively, n = display, 0 = don't:

% fig.make = [1 2 3 4 5 6 7 8 9]; % all 

% fig.make = [1 2 0 0 0 0 0 0 0]; % raw data
% fig.opt.lbl = 'raw';

% fig.make = [0 0 3 4 0 0 0 0 0]; % integrated data
% fig.opt.lbl = 'int';

fig.make = [0 0 0 0 5 0 7 8 9]; % segmentation
fig.opt.lbl = 'seg';

% ------------------------------------------------------------------------------

if diff(size(fig.make)) > 0
    fig.make = fig.make';
end

sp_order = fig.make;
i_order = find(sp_order);
sp_order(i_order) = find(i_order);

if fig.make(4)

    sp_i_offset = n_class - 1;
    sp_order = [ sp_order(1:3) ; (1:n_class)' + max(sp_order(1:3)) ; sp_order(5:end) + sp_i_offset ];
else
    
    sp_i_offset = 0;
end

n_rows = max(sp_order);
h_sp{n_rows} = [];
i_pcol = 0;

% get minima and maxima:
all_minmax = zeros(length(sp_order),2);
all_minmax(:,1) = inf;
for i_p = process_participants
    
    if fig.make(1)
        
        this_max = max(out(i_p).noint.sel_mad);
        if this_max  > all_minmax(1,2)
            
            all_minmax(1,2) = this_max + this_max/fig.opt.headroom;
        end
    end
    
    if fig.make(3)
        
        this_max = max(out(i_p).int.act_score);
        if this_max  > all_minmax(3,2)
            
            all_minmax(3,2) = this_max + this_max/fig.opt.headroom;
        end
    end
    
    if fig.make(4)
        
        for i_class = 1:n_class
            
            i_this = 3 + i_class;
            this_max = max(out(i_p).int.sel_class(:,i_class));
            if this_max  > all_minmax(i_this,2)
                
                all_minmax(i_this,2) = this_max + this_max/fig.opt.headroom;
            end
        end
    end
    
    if fig.make(9)
        
        i_this = 9 + sp_i_offset;
        this_minmax = minmax( out(i_p).int.seg_act.all(:)' );
        if this_minmax(1) < all_minmax(i_this,1)
            
            all_minmax(i_this,1) = this_minmax(1);
        end
        if this_minmax(2)  > all_minmax(i_this,2)
            
            all_minmax(i_this,2) = this_minmax(2);
        end
    end
end

if fig.opt.gray && (n_rows || isempty(info.int_param))
    
    gray_cmap.simmat = repmat( flip(0.5:0.00196:1)' ,1,3);
    gray_cmap.nov = repmat( (flip(0.9:0.000392:1)').^20,1,3);
    gray_cmap.sb = repmat( [1;0.8;0],1,3);    
    gray_cmap.seg = repmat( rescale( asin(flip(0:0.001565:1)'), 0, 1) ,1,3);
end

alphabet = 'abcdefghijklmnopqrstuvwxyz';
annotations{n_rows} = [];

monpos = get(0, 'MonitorPositions');
scrsz = monpos(fig.opt.i_mon,:);
fig_pos = scrsz;
if ~isempty(fig.opt.size_fig)
    
    fig_pos(3:4) = fig.opt.size_fig;
end
N_sp = 8 + n_class;
fig_pos(4) = round(fig_pos(4) * n_rows / N_sp );
figure('Position',fig_pos)

for i_p = process_participants
    
    i_pcol = i_pcol + 1;
    
    if fig.make(3) % Activity Score (log3 of integrated MAD + 1)
        
        plot_sel_hours_dec = out(i_p).int.sel_hours_dec;
        plot_sel_mad = out(i_p).int.act_score;

        i_sp = (sp_order(3) - 1) * n_proc_p + i_pcol;
        subplot(n_rows,n_proc_p,i_sp)
        if fig.opt.gray
            h_sp{i_sp} = plot(out(i_p).int.sel_hours_dec,plot_sel_mad,'k','LineWidth',fig.opt.ln_wdth);
        else
            h_sp{i_sp} = plot(out(i_p).int.sel_hours_dec,plot_sel_mad,'LineWidth',fig.opt.ln_wdth);
        end
        h_sp{i_sp}.Parent.Units = 'Pixels';
        h_sp{i_sp}.Parent.Position(3:4) = fig.opt.sp_size;
        box off
        line([out(i_p).int.sel_hours_dec(1),plot_sel_hours_dec(end)],[all_minmax(3,2),all_minmax(3,2)],'Color','k')
        line([out(i_p).int.sel_hours_dec(end),plot_sel_hours_dec(end)],[0,all_minmax(3,2)],'Color','k')
        ylim([0,all_minmax(3,2)])
        xlim([out(i_p).int.sel_hours_dec(1),plot_sel_hours_dec(end)])
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
        set(gca,'TickDir','out')
        ylabel('Activity Score')
        xlabel('hour')
        if fig.opt.titles
            this_title = sprintf('Activity Score for %s',info.data(i_p).descr.long);
        else
            this_title = info.data(i_p).descr.long;
        end
        title(this_title);
        if fig.opt.pnl_lbl && (i_p == process_participants(1) )
            i_row = sp_order(3);
            annotations{i_row} = annotation('textbox','String',['(',alphabet(i_row),')'],'EdgeColor','none');
            annotations{i_row}.Units = 'pixels';
            annotations{i_row}.FontSize = h_sp{i_sp}.Parent.Title.FontSize * 1.5;
            annotations{i_row}.Position = h_sp{i_sp}.Parent.Position;
            annotations{i_row}.Position(1) = annotations{i_row}.Position(1) - annotations{i_row}.Position(1)/1.5;
            annotations{i_row}.Position(4) = annotations{i_row}.Position(4)/1.5;
        end
        h_sp{i_sp}.Parent.XTickMode = 'manual';
    end
    
    if fig.make(4) % Integrated Classes
        
        for i_class = 1:n_class
            
            i_sp_order = 3 + i_class;
            i_sp = (sp_order(i_sp_order) - 1) * n_proc_p + i_pcol;
            subplot(n_rows,n_proc_p,i_sp)
            if fig.opt.gray
                h_sp{i_sp} = plot(plot_sel_hours_dec,out(i_p).int.sel_class(:,i_class),'k','LineWidth',fig.opt.ln_wdth);
            else
                h_sp{i_sp} = plot(plot_sel_hours_dec,out(i_p).int.sel_class(:,i_class),'LineWidth',fig.opt.ln_wdth);
            end
            h_sp{i_sp}.Parent.Units = 'Pixels';
            h_sp{i_sp}.Parent.Position(3:4) = fig.opt.sp_size;
            box off
            line([out(i_p).int.sel_hours_dec(1),plot_sel_hours_dec(end)],[all_minmax(i_sp_order,2),all_minmax(i_sp_order,2)],'Color','k')
            line([out(i_p).int.sel_hours_dec(end),plot_sel_hours_dec(end)],[0,all_minmax(i_sp_order,2)],'Color','k')
            ylim([0,all_minmax(i_sp_order,2)])
            xlim([out(i_p).int.sel_hours_dec(1),out(i_p).int.sel_hours_dec(end)])
            h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
            set(gca,'TickDir','out')
            ylabel('sum')
            xlabel('hour')    
            if fig.opt.titles
                this_title = sprintf('%s (%s)',info.activity_lbl{i_class},info.data(i_p).descr.long);
            else
                this_title = info.data(i_p).descr.long;
            end
            title(this_title);
            if fig.opt.pnl_lbl && (i_p == process_participants(1) )
                i_row = sp_order(i_sp_order);
                annotations{i_row} = annotation('textbox','String',['(',alphabet(i_row),')'],'EdgeColor','none');
                annotations{i_row}.Units = 'pixels';
                annotations{i_row}.FontSize = h_sp{i_sp}.Parent.Title.FontSize * 1.5;
                annotations{i_row}.Position = h_sp{i_sp}.Parent.Position;
                annotations{i_row}.Position(1) = annotations{i_row}.Position(1) - annotations{i_row}.Position(1)/1.5;
                annotations{i_row}.Position(4) = annotations{i_row}.Position(4)/1.5;
            end
            h_sp{i_sp}.Parent.XTickMode = 'manual';
        end
    end
    
    if fig.make(1) % MAD (non-integrated)
        
        plot_sel_hours_dec = info.data(i_p).hours_dec( i_adjusted );
        plot_sel_mad = out(i_p).noint.sel_mad;
        
        i_sp = (sp_order(1) - 1) * n_proc_p + i_pcol;
        subplot(n_rows,n_proc_p,i_sp) % MAD
        if fig.opt.gray
            h_sp{i_sp} = plot(plot_sel_hours_dec,plot_sel_mad,'k','LineWidth',fig.opt.ln_wdth);
        else
            h_sp{i_sp} = plot(plot_sel_hours_dec,plot_sel_mad,'LineWidth',fig.opt.ln_wdth);
        end
        h_sp{i_sp}.Parent.Units = 'Pixels';
        h_sp{i_sp}.Parent.Position(3:4) = fig.opt.sp_size;
        box off
        line([out(i_p).int.sel_hours_dec(1),plot_sel_hours_dec(end)],[all_minmax(1,2),all_minmax(1,2)],'Color','k')
        line([out(i_p).int.sel_hours_dec(end),plot_sel_hours_dec(end)],[0,all_minmax(1,2)],'Color','k')
        ylim([0,all_minmax(1,2)])
        xlim([plot_sel_hours_dec(1),plot_sel_hours_dec(end)])
        hour_step = diff(h_sp{i_sp}.Parent.XTick(1:2));
        h_sp{i_sp}.Parent.XTick = old_XTick;        
        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
        set(gca,'TickDir','out')
        ylabel('MAD')
        xlabel('hour')
        if fig.opt.titles
            this_title = sprintf('MAD for %s',info.data(i_p).descr.long);
        else
            this_title = info.data(i_p).descr.long;
        end
        title(this_title);
        if fig.opt.pnl_lbl && (i_p == process_participants(1) )
            i_row = sp_order(1);
            annotations{i_row} = annotation('textbox','String',['(',alphabet(i_row),')'],'EdgeColor','none');
            annotations{i_row}.Units = 'pixels';
            annotations{i_row}.FontSize = h_sp{i_sp}.Parent.Title.FontSize * 1.5;
            annotations{i_row}.Position = h_sp{i_sp}.Parent.Position;
            annotations{i_row}.Position(1) = annotations{i_row}.Position(1) - annotations{i_row}.Position(1)/1.5;
            annotations{i_row}.Position(4) = annotations{i_row}.Position(4)/1.5;
        end
        h_sp{i_sp}.Parent.XTickMode = 'manual';
    end
    
    if fig.make(2) % Activity States (classes, non-integrated)
        
        i_sp = (sp_order(2) - 1) * n_proc_p + i_pcol;
        subplot(n_rows,n_proc_p,i_sp) 
        disp_sel_class = out(i_p).noint.sel_class;
        disp_sel_class( disp_sel_class == 0) = NaN;
        if fig.opt.gray
            h_sp{i_sp} = plot(plot_sel_hours_dec,disp_sel_class,'k.','MarkerSize',fig.opt.mk_size);
        else
            h_sp{i_sp} = plot(plot_sel_hours_dec,disp_sel_class,'.','MarkerSize',fig.opt.mk_size);
        end
        h_sp{i_sp}.Parent.Units = 'Pixels';
        h_sp{i_sp}.Parent.Position(3:4) = fig.opt.sp_size;
        box off
        line([out(i_p).int.sel_hours_dec(1),plot_sel_hours_dec(end)],[max(out(i_p).noint.sel_class)+0.5,max(out(i_p).noint.sel_class)+0.5],'Color','k')
        line([out(i_p).int.sel_hours_dec(end),plot_sel_hours_dec(end)],[0.5, max(out(i_p).noint.sel_class)+0.5],'Color','k')
        xlim([plot_sel_hours_dec(1),plot_sel_hours_dec(end)])
        ylim([0.5, max(out(i_p).noint.sel_class)+0.5]) % exclude zeros (formerly NaN) from axis
        h_sp{i_sp}.Parent.XTick = old_XTick;
        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
        these_class = unique(out(i_p).noint.sel_class);
        these_class( these_class == 0) = [];
        h_sp{i_sp}.Parent.YTick = these_class;
        h_sp{i_sp}.Parent.YTickLabel = info.activity_lbl(these_class);
        set(gca,'TickDir','out')
        xlabel('hour')
        if fig.opt.titles
            this_title = sprintf('Activity States for %s',info.data(i_p).descr.long);
        else
            this_title = info.data(i_p).descr.long;
        end
        title(this_title);
        if fig.opt.pnl_lbl && (i_p == process_participants(1) )
            i_row = sp_order(2);
            annotations{i_row} = annotation('textbox','String',['(',alphabet(i_row),')'],'EdgeColor','none');
            annotations{i_row}.Units = 'pixels';
            annotations{i_row}.FontSize = h_sp{i_sp}.Parent.Title.FontSize * 1.5;
            annotations{i_row}.Position = h_sp{i_sp}.Parent.Position;
            annotations{i_row}.Position(1) = annotations{i_row}.Position(1) - annotations{i_row}.Position(1)/1.5;
            annotations{i_row}.Position(4) = annotations{i_row}.Position(4)/1.5;
        end
        h_sp{i_sp}.Parent.XTickMode = 'manual';
    end

    if fig.make(5) % Activity Self-Similarity 
               
        i_sp = (sp_order(5 + sp_i_offset) - 1) * n_proc_p + i_pcol;
        subplot(n_rows,n_proc_p,i_sp)
        h_sp{i_sp} = imagesc(out(i_p).int.sim_mat);
        h_sp{i_sp}.Parent.Units = 'Pixels';
        h_sp{i_sp}.Parent.Position(3:4) = fig.opt.sp_size;
        box off
        these_XLim = get(gca,'XLim');
        these_YLim = get(gca,'YLim');
        line(these_XLim,[these_YLim(1),these_YLim(1)],'Color','k')
        line([these_XLim(2),these_XLim(2)],these_YLim,'Color','k')
        h_sp{i_sp}.Parent.YTick = [];
        h_sp{i_sp}.Parent.XTick = i_map_xticks;
        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
        set(gca,'TickDir','out')
        xlabel('hour')
        if fig.opt.titles
            this_title = sprintf('Activity Self-Similarity for %s',info.data(i_p).descr.long);
        else
            this_title = info.data(i_p).descr.long;
        end
        title(this_title);
        if fig.opt.gray
            colormap(h_sp{i_sp}.Parent,gray_cmap.simmat)
        end
        if fig.opt.pnl_lbl && (i_p == process_participants(1) )
            i_row = sp_order(5 + sp_i_offset);
            annotations{i_row} = annotation('textbox','String',['(',alphabet(i_row),')'],'EdgeColor','none');
            annotations{i_row}.Units = 'pixels';
            annotations{i_row}.FontSize = h_sp{i_sp}.Parent.Title.FontSize * 1.5;
            annotations{i_row}.Position = h_sp{i_sp}.Parent.Position;
            annotations{i_row}.Position(1) = annotations{i_row}.Position(1) - annotations{i_row}.Position(1)/1.5;
            annotations{i_row}.Position(4) = annotations{i_row}.Position(4)/1.5;
        end
        h_sp{i_sp}.Parent.XTickMode = 'manual';
    end
    
    if fig.make(6) % Activity Novelty 
        
        i_sp = (sp_order(6 + sp_i_offset) - 1) * n_proc_p + i_pcol;
        subplot(n_rows,n_proc_p,i_sp)
        h_sp{i_sp} = imagesc(out(i_p).int.nov);
        h_sp{i_sp}.Parent.Units = 'Pixels';
        h_sp{i_sp}.Parent.Position(3:4) = fig.opt.sp_size;
        box off
        these_XLim = get(gca,'XLim');
        these_YLim = get(gca,'YLim');
        line(these_XLim,[these_YLim(1),these_YLim(1)],'Color','k')
        line([these_XLim(2),these_XLim(2)],these_YLim,'Color','k')
        h_sp{i_sp}.Parent.XTick = i_map_xticks;
        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
        set(gca,'TickDir','out')
        ylabel('granularity')
        xlabel('hour')
        if fig.opt.titles
            this_title = sprintf('Activity Novelty for %s',info.data(i_p).descr.long);
        else
            this_title = info.data(i_p).descr.long;
        end
        title(this_title);
        if fig.opt.gray
            colormap(h_sp{i_sp}.Parent,gray_cmap.nov)
        end
        if fig.opt.pnl_lbl && (i_p == process_participants(1) )
            i_row = sp_order(6 + sp_i_offset);
            annotations{i_row} = annotation('textbox','String',['(',alphabet(i_row),')'],'EdgeColor','none');
            annotations{i_row}.Units = 'pixels';
            annotations{i_row}.FontSize = h_sp{i_sp}.Parent.Title.FontSize * 1.5;
            annotations{i_row}.Position = h_sp{i_sp}.Parent.Position;
            annotations{i_row}.Position(1) = annotations{i_row}.Position(1) - annotations{i_row}.Position(1)/1.5;
            annotations{i_row}.Position(4) = annotations{i_row}.Position(4)/1.5;
        end
        h_sp{i_sp}.Parent.XTickMode = 'manual';
    end
    
    if fig.make(7) % Activity Boundaries (novelty peaks)
        
        i_sp = (sp_order(7 + sp_i_offset) - 1) * n_proc_p + i_pcol;
        subplot(n_rows,n_proc_p,i_sp)
        h_sp{i_sp} = imagesc(out(i_p).int.i_seg.all);
        h_sp{i_sp}.Parent.Units = 'Pixels';
        h_sp{i_sp}.Parent.Position(3:4) = fig.opt.sp_size;
        box off
        these_XLim = get(gca,'XLim');
        these_YLim = get(gca,'YLim');
        line(these_XLim,[these_YLim(1),these_YLim(1)],'Color','k')
        line([these_XLim(2),these_XLim(2)],these_YLim,'Color','k')
        h_sp{i_sp}.Parent.XTick = i_map_xticks;
        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
        set(gca,'TickDir','out')
        ylabel('granularity')
        xlabel('hour')
        if fig.opt.titles
            this_title = sprintf('Activity Boundaries for %s',info.data(i_p).descr.long);
        else
            this_title = info.data(i_p).descr.long;
        end
        title(this_title);
        if fig.opt.gray
            colormap(h_sp{i_sp}.Parent,gray_cmap.sb)
        end
        if fig.opt.pnl_lbl && (i_p == process_participants(1) )
            i_row = sp_order(7 + sp_i_offset);
            annotations{i_row} = annotation('textbox','String',['(',alphabet(i_row),')'],'EdgeColor','none');
            annotations{i_row}.Units = 'pixels';
            annotations{i_row}.FontSize = h_sp{i_sp}.Parent.Title.FontSize * 1.5;
            annotations{i_row}.Position = h_sp{i_sp}.Parent.Position;
            annotations{i_row}.Position(1) = annotations{i_row}.Position(1) - annotations{i_row}.Position(1)/1.5;
            annotations{i_row}.Position(4) = annotations{i_row}.Position(4)/1.5;
        end
        h_sp{i_sp}.Parent.XTickMode = 'manual';
    end
    
    if fig.make(8) % Rectified and Reduced Activity Boundaries
        
        i_sp = (sp_order(8 + sp_i_offset) - 1) * n_proc_p + i_pcol;
        subplot(n_rows,n_proc_p,i_sp)
        h_sp{i_sp} = imagesc( out(i_p).int.i_seg.all_rect_rdx );
        h_sp{i_sp}.Parent.Units = 'Pixels';
        h_sp{i_sp}.Parent.Position(3:4) = fig.opt.sp_size;
        box off
        these_XLim = get(gca,'XLim');
        these_YLim = get(gca,'YLim');
        line(these_XLim,[these_YLim(1),these_YLim(1)],'Color','k')
        line([these_XLim(2),these_XLim(2)],these_YLim,'Color','k')
        h_sp{i_sp}.Parent.XTick = i_map_xticks;
        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
        set(gca,'TickDir','out')
        ylabel('granularity')
        xlabel('hour')
        if fig.opt.titles
            this_title = sprintf('Rectified and Reduced Activity Boundaries for %s',info.data(i_p).descr.long);
        else
            this_title = info.data(i_p).descr.long;
        end
        title(this_title);
        if fig.opt.gray
            colormap(h_sp{i_sp}.Parent,gray_cmap.sb)
        end
        if fig.opt.pnl_lbl && (i_p == process_participants(1) )
            i_row = sp_order(8 + sp_i_offset);
            annotations{i_row} = annotation('textbox','String',['(',alphabet(i_row),')'],'EdgeColor','none');
            annotations{i_row}.Units = 'pixels';
            annotations{i_row}.FontSize = h_sp{i_sp}.Parent.Title.FontSize * 1.5;
            annotations{i_row}.Position = h_sp{i_sp}.Parent.Position;
            annotations{i_row}.Position(1) = annotations{i_row}.Position(1) - annotations{i_row}.Position(1)/1.5;
            annotations{i_row}.Position(4) = annotations{i_row}.Position(4)/1.5;
        end
        h_sp{i_sp}.Parent.XTickMode = 'manual';
    end
    
    if fig.make(9) % Segmented Activity Score
        
        i_sp_order = 9 + sp_i_offset;
        i_sp = (sp_order(i_sp_order) - 1) * n_proc_p + i_pcol;
        subplot(n_rows,n_proc_p,i_sp)
        h_sp{i_sp} = imagesc( out(i_p).int.seg_act.all );
        h_sp{i_sp}.Parent.Units = 'Pixels';
        h_sp{i_sp}.Parent.Position(3:4) = fig.opt.sp_size;
        box off
        these_XLim = get(gca,'XLim');
        these_YLim = get(gca,'YLim');
        line(these_XLim,[these_YLim(1),these_YLim(1)],'Color','k')
        line([these_XLim(2),these_XLim(2)],these_YLim,'Color','k')
        caxis([ all_minmax(i_sp_order,1) ,all_minmax(i_sp_order,2)])
        h_sp{i_sp}.Parent.XTick = i_map_xticks;
        h_sp{i_sp}.Parent.XTickLabel = new_XTick_lbl;
        set(gca,'TickDir','out')
        ylabel('granularity')
        xlabel('hour')
        if fig.opt.titles
            this_title = sprintf('Segmented Activity Score for %s',info.data(i_p).descr.long);
        else
            this_title = info.data(i_p).descr.long;
        end
        title(this_title);
        if fig.opt.gray
            colormap(h_sp{i_sp}.Parent,gray_cmap.seg)
        end
        if fig.opt.pnl_lbl && (i_p == process_participants(1) )
            i_row = sp_order(i_sp_order);
            annotations{i_row} = annotation('textbox','String',['(',alphabet(i_row),')'],'EdgeColor','none');
            annotations{i_row}.Units = 'pixels';
            annotations{i_row}.FontSize = h_sp{i_sp}.Parent.Title.FontSize * 1.5;
            annotations{i_row}.Position = h_sp{i_sp}.Parent.Position;
            annotations{i_row}.Position(1) = annotations{i_row}.Position(1) - annotations{i_row}.Position(1)/1.5;
            annotations{i_row}.Position(4) = annotations{i_row}.Position(4)/1.5;
        end
        h_sp{i_sp}.Parent.XTickMode = 'manual';
    end
end

if ~isempty(fig.opt.lbl)
    
    fig_ext = 'tif';
    full_fn_fig = [fig.opt.savepath,'/FIG_',fig.opt.lbl,'_',info.dh_lbl,'.',fig_ext];
    set(gcf,'paperunits','normalized');
    paper_pos = get(gcf,'paperposition');
    paper_pos(3) = paper_pos(3) * fig.opt.h_resc; 
    set(gcf,'paperposition',paper_pos);
    print(gcf,full_fn_fig,'-dtiffn',sprintf('-r%i',fig.opt.dpi));
end

disp(' ')
disp('done')