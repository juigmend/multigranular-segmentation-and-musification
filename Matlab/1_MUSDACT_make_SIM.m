%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                                 SIMULATION OF                                %
%                            MEAN ABSOLUTE DEVIATION                           %
%                       OF THIGH-WORN ACCELEROMETER DATA                       %
%                                                                              %
%                               25 OCTOBER, 2022                               %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tested with Matlab R2017b

% Instructions:
%     Edit the values indicated with an arrow like this: <--- (length of the arrow may vary)
%     Run the program, close your eyes and hope for the best.

% ==============================================================================

   paths_data  = '';  % <--- path for data (raw and processed)

 classes_names = {'walking','upright','sitting','lying'}; % <---
 class_weights = [1,1,-0.5,0]; % <---
  
     epoch_p_h = 60 * 12; % 5s. epoch
     
gausswin_alpha = 4; % walking, sitting
    tukeywin_r = 0.5; % upright
    
% ------------------------------------------------------------------------------

win_length = epoch_p_h;

for i_p = 1:2
    
    % ..............................................................................
    
    if i_p == 1
        
        participant_description = 'LOW';
        classif_mtx(1,:) = [ 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 1 ]; % <--- walking
        classif_mtx(2,:) = [ 0 0 0 0 0 0 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ]; % <--- upright
        classif_mtx(3,:) = [ 0 1 0 1 0 0 1 1 1 1 1 0 0 1 1 1 1 0 1 0 0 1 0 1 ]; % <--- sitting
        classif_mtx(4,:) = [ 1 1 1 1 1 1 0 1 0 0 1 1 1 0 0 0 1 1 1 1 1 1 1 1 ]; % <--- lying
        % FOO_time_index = [ 4 5 6 7 8 9 0 1 2 1 2 3 4 5 6 7 8 9 0 1 2 1 2 3 ]; % hours
        
    elseif i_p == 2
        
        participant_description = 'HIGH';
        classif_mtx(1,:) = [ 1 0 0 1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 ]; % <--- walking
        classif_mtx(2,:) = [ 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 ]; % <--- upright
        classif_mtx(3,:) = [ 1 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 0 0 1 0 1 ]; % <--- sitting
        classif_mtx(4,:) = [ 1 1 1 1 1 1 0 1 0 0 0 0 0 1 0 0 1 1 1 1 1 1 1 1 ]; % <--- lying
        % FOO_time_index = [ 4 5 6 7 8 9 0 1 2 1 2 3 4 5 6 7 8 9 0 1 2 1 2 3 ]; % hours
        
    end
  
    % ------------------------------------------------------------------------------
  
    now_unix = int64( posixtime( datetime('now') ));
    n_classes = length(classes_names);
    
    if n_classes ~= size(classif_mtx,1)
        error('n_classes ~= size(classif_mtx,1)')
    end
    
    time_length = size(classif_mtx,2);
    
    total_time_length =  time_length * epoch_p_h;
    
    activities.tStamps = ( now_unix : 5000 : now_unix + ( total_time_length * 5000 ) - 5000 )'; % 5s. epoch
    
    classif_mtx_w = classif_mtx;
    classif_bool = logical(classif_mtx);
    classif_vec = zeros(1,time_length);
    classif_vec_w = classif_vec;
    token_class = 1;
    
    for i_class = flip(1:n_classes)
        
        i_this = classif_bool(i_class,:);
        
        classif_vec(i_this)  = classif_bool(i_class,i_this) * token_class;
        token_class = token_class + 1;
        
        classif_mtx_w(i_class,i_this) = classif_bool(i_class,i_this) * class_weights(i_class);
    end
    
    classif_vec_NaN = classif_vec;
    classif_vec_NaN(classif_vec_NaN == 0) = NaN;
    classif_vec_NaN = classif_vec_NaN - 1;
    classif_vec_NaN_expanded = zeros(1,total_time_length);
    this_start = 1;
    
    for i_time = 1:time_length
        
        this_end = this_start + win_length - 1;
        classif_vec_NaN_expanded(this_start:this_end) = classif_vec_NaN(i_time);
        this_start = this_end+1;
    end
    
    activities.classification = classif_vec_NaN_expanded';
    
    classif_mtx_w( classif_mtx_w < 0) = 0;
    classif_mtx_w_expanded = zeros(n_classes,total_time_length);
    win_length_half = floor(win_length/2);
    
    for i_class = 1:n_classes
        
        this_start = 1;
        
        for i_time = 1:time_length
            
            this_end = this_start + win_length - 1;
            
            classif_mtx_w_expanded(i_class,this_start+win_length_half) = classif_mtx_w(i_class,i_time);
            
            if (i_class == 1) || (i_class == 3) % 'walking' , 'sitting'
                
                classif_mtx_w_expanded(i_class,this_start:this_end) = ...
                    conv(...
                    classif_mtx_w_expanded(i_class,this_start:this_end) , ...
                    gausswin(win_length,gausswin_alpha) ...
                    , 'same');
                
            elseif i_class == 2 % 'upright'
                
                classif_mtx_w_expanded(i_class,this_start:this_end) = ...
                    conv(...
                    classif_mtx_w_expanded(i_class,this_start:this_end) , ...
                    tukeywin(win_length,tukeywin_r)  ...
                    , 'same');
                
            end
            
            % no 'lying' :-P
            
            this_start = this_end+1;
        end
    end
    
    MAD = sum(classif_mtx_w_expanded)';
    MAD = MAD / max(MAD);
    MAD = MAD * 0.7; % rescale
    activities.ukk.mad = MAD;
    
    this_fname = [paths_data,'/SIM_',participant_description,'.mat'];
    
    save(this_fname,'activities')
end
