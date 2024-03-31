function [data, expmnt] = ssvepCFS
% Info
% 
% Pending:
% [X] Custom adjustment of x_displacement.
% [X] Present masks depending on which frame we are.
% [X] Present stimuli depending on which frame we are.
% [X] Add fade-in & fade-out for each stimulus.
% [X] Test keyboard check (do we have problems with timing?).
%     Apparently there is no lag!!
% [X] Save to Github.
% [X] Add overall fade-in & fade-out for entire block.
% [ ] Verify frames presented: save conditions vectors.
% [ ] Adjust mask and faces sizes: based on visual angle degrees.
% [ ] Change place of Hz estimations (after PsychToolbox start).
% [ ] Save mask, target, vergence bars positions: expmnt.

clc;
fprintf(    '_________________________________________________________________________\n');
fprintf('\n\n      CFS & ssVEP experiment   \n\n');
fprintf(    '      Code: Manuel Mejia   \n');
fprintf(    '      Date: 2024-March-29   \n');
fprintf('\n\n_________________________________________________________________________\n');

expmnt.savescreen = 0;
addpath('C:\Users\cthul\Dropbox\PROYECTOS\0-PhD\MLtoolboxes\SHINEtoolbox');

%% 60 or 85 Hz
% Frequencies that can be tagged depending on monitor screen rate:
% 1./([1:50].*(1000./60)).*1000 % 1.2, 6, 8.5714 = 50, 10, 7
% 1./([50, 10, 7].*(1000./60)).*1000
% If we define Hz for masks and baseline stimuli:
% 1./([1:80].*(1000./85)).*1000 % 1.2143, 6.0714, 8.5  = 70, 14, 10
% 1./([70, 14, 10].*(1000./85)).*1000

%% Set experiment variables
dir_images = '../stimuli/caras_experimento/';

expmnt.masks_hz =       8.5; % approx, is defined according to monitor framerate
expmnt.baseline_hz =    6;
expmnt.trial_duration = 60;
expmnt.fade_in_time =   2;
expmnt.fade_out_time =  2;
expmnt.odd_frequency =  5; % each 5th baseline stim, the oddball

expmnt.monitor_hz = 60;
mask_framerate = round(expmnt.monitor_hz/expmnt.masks_hz);
baseline_framerate = round(expmnt.monitor_hz/expmnt.baseline_hz);
oddball_framerate = baseline_framerate*expmnt.odd_frequency;

expmnt.stimuli_width = 150;
expmnt.stimuli_height = 300;

% Masks box:
expmnt.mask_width = 400;
expmnt.mask_height = 400;
expmnt.num_masks = 500;

expmnt.im_white    = 255;
expmnt.img_contrast_sd         = 12; % approx 20% Michelson contrast. Need to check. %%%  EDIT BEFORE EXP. %%% 
expmnt.cfs_masks_contrast_sd   = 82; % the mean of the original function make_mondrian_masks
expmnt.gaussian_envelope   = 0; % If blur edges with gaussian envelope.

expmnt.stimuli_width       = 100; %150x200; Best: 120x160, for CFS smller than mask, and for having enough screen for crowding.
expmnt.stimuli_height      = 133;

expmnt.position_task_offset    = 0;
expmnt.vergence_bar = [0, 0, 32, 240];
expmnt.show_vergence_bars = 1;
expmnt.vergence_bar_file       = '../stimuli/texture_vertical.bmp';
expmnt.dir_of_stereograms     = '../stimuli/stereograms';
expmnt.stereogram_shapes      = {'square', 'T', 'square', 'square', 'upper-square', 'lower-square', 'two-squares', 'L', 'U', 'big-square'};
expmnt.stereogram_for_pf               = [2, 8, 9, 2, 8];
expmnt.stereogram_for_exp              = [2, 8, 9, 2, 8];

expmnt.fix_color = [1,1,1];
expmnt.img_lum_mean    = 120; % background x im_white, equal to background mean luminance.
expmnt.background = 0.47;
expmnt.cfs_with_one_mirror_set = 0; % Set to 1 if only using mirrors for one eye: for using Eyelink tracking.
expmnt.x_displacement = 300;
if expmnt.cfs_with_one_mirror_set == 1
    expmnt.cfs_1mirror_set_left_prop_of_displace       = 0.90;
    expmnt.cfs_1mirror_set_right_prop_of_displace      = 1.10;
    expmnt.cfs_1mirror_set_size_ratio                  = 0.86;
else
    expmnt.cfs_1mirror_set_left_prop_of_displace       = 1.00;
    expmnt.cfs_1mirror_set_right_prop_of_displace      = 1.00;
    expmnt.cfs_1mirror_set_size_ratio                  = 1.00;
end

expmnt.PTB_SkipSyncTexts_config = 1;
expmnt.debug_pt = 0;
expmnt.BeampositionQueryWorkaround = 1;
expmnt.set_custom_screen_size = 0;



%% Estimated vars within experiment:

target_box = [0 0 expmnt.stimuli_width expmnt.stimuli_height]; % Position of target.
% Masks box:
mask_box = [0, 0, expmnt.mask_width, expmnt.mask_height];
masks_time_interval = 1/expmnt.masks_hz;
% Make the contingency window, to check if eyegaze is here:
gazeRect = [0 0 200 200];
% Fixation cross specs:
fix_color = expmnt.fix_color;
fix_size = 10;
% Line width for fixation cross:
fix_lineWidth = 2;

font_size = 15;

which_stereogram    = expmnt.stereogram_for_exp;

pos_first_text_line = 0.40;

%% Prepare stimuli

% Get list of filenames of stimuli:
familiar_images = dir([dir_images filesep 'familiar*.png']);
familiar_images = {familiar_images.name};
unfamiliar_images = dir([dir_images filesep 'unfamiliar*.png']);
unfamiliar_images = {unfamiliar_images.name};
assert(length(familiar_images)>0 & length(unfamiliar_images)>0);

%% Prepare conditions vectors

% Frames for an entire trial:
frames_conds = ones(1, expmnt.trial_duration*expmnt.monitor_hz);

% Faces vector:
frames_faces = nan(1, length(frames_conds)/baseline_framerate);
temp = [];
for i = 1:ceil(length(frames_faces)/length(unfamiliar_images))
    temp = [temp, Shuffle(1:length(unfamiliar_images))];
end
% This is a vector with indices for each face to be presented:
frames_faces = temp(1:length(frames_faces));
clearvars temp;
% Every Xth image, present the oddball stimulus:
idx_oddball = expmnt.odd_frequency:expmnt.odd_frequency:length(frames_faces);
frames_faces(idx_oddball) = Shuffle((1:length(familiar_images))+length(unfamiliar_images));

% Expand each image index to the number of frames to be presented:
frames_faces = Expand(frames_faces, baseline_framerate, 1);




%% Start PsychToolbox

% Call some default settings:
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', expmnt.PTB_SkipSyncTexts_config); % This is to avoid 'errors' when VBL does not work correctly.
if expmnt.debug_pt
        PsychDebugWindowConfiguration(0, 0.5); %0.4 is 40% opaque.
end
if expmnt.BeampositionQueryWorkaround
    Screen('Preference', 'ConserveVRAM', 4096); % Setting kPsychUseBeampositionQueryWorkaround
end
% Check how many screens are there (0 means 1):
screens = Screen('Screens');
% Draw the external screen if available:
screenNumber = max(screens);
maxPriorityLevel = MaxPriority(screenNumber); 
if expmnt.set_custom_screen_size
    oldResolution = Screen('Resolution', screenNumber, expmnt.screenX, expmnt.screenY);
end
% Define black and white indices:
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
% 	gray = (white + black)/2;
% Set coordinates for fixation:
fix_xCoords = [-fix_size, fix_size, 0,           0];
fix_yCoords = [0,         0,        -fix_size,   fix_size];
fixCoords = [fix_xCoords; fix_yCoords];
fixCoordsLeft = fixCoords*expmnt.cfs_1mirror_set_size_ratio;
fixCoordsRight = fixCoords;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set available keys to press:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get key codes:
escapeKey =     KbName('ESCAPE');
spaceKey =      KbName('space');
key1 =          KbName('LeftArrow');
key2 =          KbName('RightArrow');
responseKeys =  zeros(1, 256);
responseKeys([escapeKey, spaceKey, key1, key2]) = 1;
pauseKeys =     zeros(1, 256);
pauseKeys(spaceKey) = 1;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start PsychToolbox window:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open an on screen window:
[w, wRect] = PsychImaging('OpenWindow', screenNumber, expmnt.background);
% Get size of on screen window:
[screenXpixels, screenYpixels] = Screen('WindowSize', w);
% Query the frame duration. It will be used to control the while loop.
ifi = Screen('GetFlipInterval', w);
% Get the centre coordinate of the window.
[xCenter, yCenter] = RectCenter(wRect);
% Set up alpha-blending for smooth (anti-aliased) lines????
Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('HideCursorHelper', w);


%% Duty cycle for sinusoidal modulation
Hz = 60;
flickerRate = expmnt.baseline_hz; % 6 Hz
cycle = round(linspace(0,360,Hz));                  %round(1/(frames)*1000)));
cycle = repmat(cycle,[1,expmnt.trial_duration,1]);
trans = cosd(cycle*flickerRate+180)*0.5+0.5;        %Create sinusoidal levels for pixel level
% trans(end) = [];                                    %Remove last element (as first is the same)
%180 = shift phase; *0.5+0.5 = Transparency level 0-1
%mins = findminima(trans);
flipped = 1 - trans;
[peaks,locs]=findpeaks(flipped);      %Flip then find troughs
mins = locs;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vergence bar:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vergence(1).img = imread(expmnt.vergence_bar_file);
vergence(1).img = mat2gray(vergence(1).img, [0, 3]);
if expmnt.show_vergence_bars == 0
    vergence(1).img(:, :, 2) = zeros(size(vergence(1).img, 1), size(vergence(1).img, 2));
end
vergence(1).tex = Screen('MakeTexture', w, vergence(1).img);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load stereograms:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stereograms_list = dir(fullfile(expmnt.dir_of_stereograms, '*l.jpg'));
for i = 1:length(stereograms_list)
    stereograms(i).img_left =  imread(fullfile(expmnt.dir_of_stereograms, sprintf('stereogram%02d_l.jpg', i)));
    stereograms(i).img_right = imread(fullfile(expmnt.dir_of_stereograms, sprintf('stereogram%02d_r.jpg', i)));

    stereograms(i).tex_left =  Screen('MakeTexture', w, stereograms(i).img_left);
    stereograms(i).tex_right = Screen('MakeTexture', w, stereograms(i).img_right);
end
clearvars stereograms_list;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimuli fixed positions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that generates positions for images that remain static throughout
% entire experiment: 
% 1) vergence bars, 2) Fixation cross, 3) textboxes, 4) stereograms
[vergence_bars_positions, fix_position_left, fix_position_right, textBox, stereogram_box] =...
    create_static_images_boxes(expmnt, xCenter, yCenter, screenYpixels, mask_box);

% Center ContingencyWindow:
gazeRect = CenterRectOnPointd(gazeRect, xCenter, yCenter);

% Set size of masks & target depending if resize needed on left eye:
mask_box_left =     mask_box * expmnt.cfs_1mirror_set_size_ratio;
mask_box_right =    mask_box;
target_box_left =   target_box * expmnt.cfs_1mirror_set_size_ratio;
target_box_right =  target_box;

% Set masks & target to the center of the screen:
mask_box_left =     CenterRectOnPointd(mask_box_left, xCenter, yCenter);
mask_box_right =    CenterRectOnPointd(mask_box_right, xCenter, yCenter);
target_box_left =   CenterRectOnPointd(target_box_left, xCenter, yCenter);
target_box_right =  CenterRectOnPointd(target_box_right, xCenter, yCenter);

% When task position needed, set amount of offset for target stimuli:
left_pos_offset =   [-expmnt.position_task_offset, 0, -expmnt.position_task_offset, 0];
right_pos_offset =  [ expmnt.position_task_offset, 0,  expmnt.position_task_offset, 0];



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First screen: loading images. Text: 'Preparando las imagenes...'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Screen('TextSize', w, font_size);
KbQueueCreate([], responseKeys);
KbQueueStart;
KbQueueFlush();
x_displacement_set = 0;
while x_displacement_set == 0
    Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
    Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
    % DrawFormattedText(w, 'Preparando las imagenes...', 'center', screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.left);
    DrawFormattedText(w, 'Preparando las imagenes...', 'center', screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.right);
    Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
    Screen('Flip', w);
    [pressed, firstPress] = KbQueueCheck;
    if pressed
        if firstPress(key1)
            expmnt.x_displacement = expmnt.x_displacement - 5;
        end
        if firstPress(key2)
            expmnt.x_displacement = expmnt.x_displacement + 5;
        end
        if firstPress(spaceKey)
            x_displacement_set = 1;
        end
        [vergence_bars_positions, fix_position_left, fix_position_right, textBox, stereogram_box] =...
            create_static_images_boxes(expmnt, xCenter, yCenter, screenYpixels, mask_box);
    end
end
if expmnt.savescreen; printscreenArray{1} = Screen('GetImage', w); end; %% Printscreen and save file.




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and edit stimuli for this participant:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit images: lumMatch, background colour, and convert into textures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All faces
if 1
    [loadedFacesEdited, loadedFaces, mask_frgd, mask_bkgd] =...
        load_edit_convert_imgs([unfamiliar_images, familiar_images], dir_images, expmnt, w);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create and prepare masks:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_mondrian_masks(sz_x, sz_y, n_masks, shape, selection)
% shape=2, circles. % selection=2, grayscale.
expmnt.num_masks = ceil(length(frames_conds)/mask_framerate); % Minimum
temp = Expand(1:expmnt.num_masks, mask_framerate, 1);
frames_masks = temp(1:length(frames_conds));
clearvars temp;
masks_raw = make_mondrian_masks(expmnt.mask_width, expmnt.mask_height, expmnt.num_masks, 2, 2);
% Same luminance of masks and targets. Control contrast:
masks = lumMatch(masks_raw, [], [expmnt.img_lum_mean, expmnt.cfs_masks_contrast_sd]); 
% Convert masks into textures:
for i = 1:length(masks)
    % Make masks into textures:
    masksTex(i).all = Screen('MakeTexture', w, masks{i});
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second screen: images loaded. Text: 'Imagenes preparadas.'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DrawFormattedText(w, 'Imagenes listas.', 'center', screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.right);
Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
Screen('Flip', w);
if expmnt.savescreen; printscreenArray{end+1} = Screen('GetImage', w); end; %% Printscreen and save file.
WaitSecs(.1);
Screen('HideCursorHelper', w);
if expmnt.debug_pt == 1
    WaitSecs(2);
else
    KbWait;
    WaitSecs(.1);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third Screen: stereogram to confirm vergence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw stereograms:
Screen('DrawTextures', w, stereograms(which_stereogram(1)).tex_left, [], stereogram_box.left, []);
Screen('DrawTextures', w, stereograms(which_stereogram(1)).tex_right, [], stereogram_box.right, []);
Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
Screen('Flip', w);
WaitSecs(.1);
if expmnt.savescreen; printscreenArray{end+1} = Screen('GetImage', w); end; %% Printscreen and save file.
if expmnt.debug_pt == 1
    WaitSecs(2);
else
    KbWait;
    WaitSecs(.1);
end


%% _________________________________________________________________________
%|
%|
%|
%|
%|                    BEGINNING OF EXPERIMENT 
%|
%|
%|
%|__________________________________________________________________________


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% screen: press any key to begin. 'Presiona cualquier tecla \n\npara iniciar.'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DrawFormattedText(w, 'Presione cualquier tecla para iniciar.', 'center',...
    screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.right);
Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
Screen('Flip', w);
if expmnt.savescreen; printscreenArray{end+1} = Screen('GetImage', w); end; %% Printscreen and save file.
if expmnt.debug_pt == 1
    WaitSecs(2);
else
    KbWait;
    WaitSecs(.5);
end


%% Block loop
block_end = 0;
timestamps = nan(1, length(frames_conds));
trial_responses = [];
trial_response_times = [];
trial_response = [];
trial_response_time = [];
KbQueueFlush();

% Mask position and potency:
this_trial__xy_target_masks_displacement = ...
            [expmnt.x_displacement, 0, expmnt.x_displacement, 0];
this_trial__masks_position_xy =...
    mask_box_right + expmnt.cfs_1mirror_set_right_prop_of_displace * this_trial__xy_target_masks_displacement;
this_trial__mask_alpha = 1;

% Target position:
this_trial__target_position_xy = ...
    target_box_left - expmnt.cfs_1mirror_set_right_prop_of_displace * this_trial__xy_target_masks_displacement;

% Fade in:
fade_in_frames = expmnt.fade_in_time * expmnt.monitor_hz;
fade_in_cycle = round(linspace(0, this_trial__mask_alpha, fade_in_frames), 2);
temp = Shuffle(1:expmnt.num_masks);
temp = Expand(temp, mask_framerate, 1);
fade_in_masks = temp(1:length(fade_in_cycle));
clearvars temp;
for frame = 1:length(fade_in_cycle)
    Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
    % Draw masks:
    Screen('DrawTextures', w, masksTex(fade_in_masks(frame)).all, [], this_trial__masks_position_xy, 180, [], fade_in_cycle(frame));
    Screen('Flip', w);
end

for frame = 1:length(frames_conds)
    
    Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
    % Draw target:
    Screen('DrawTextures', w, loadedFaces(frames_faces(frame)).theTexture, [], this_trial__target_position_xy, 0, [], 1-flipped(frame));
    % Draw masks:
    Screen('DrawTextures', w, masksTex(frames_masks(frame)).all, [], this_trial__masks_position_xy, 180, [], this_trial__mask_alpha);
    
    % Gather response:
    [pressed, firstPress] = KbQueueCheck;
    if pressed
        if firstPress(key1)
            trial_response = 0;
            trial_response_time = firstPress(key1);
        end
        if firstPress(key2)
            trial_response = 1;
            trial_response_time = firstPress(key2);
        end
        if firstPress(spaceKey)
            block_end = 1;
            break;
        end
        trial_responses = [trial_responses, trial_response];
        trial_response_times = [trial_response_times, trial_response_time];
    end
    
    vbl = Screen('Flip', w);
    timestamps(frame) = vbl;
end

% Fade out:
fade_out_frames = expmnt.fade_out_time * expmnt.monitor_hz;
fade_out_cycle = round(linspace(this_trial__mask_alpha, 0, fade_out_frames), 2);
temp = Shuffle(1:expmnt.num_masks);
temp = Expand(temp, mask_framerate, 1);
fade_out_masks = temp(1:length(fade_out_cycle));
for frame = 1:length(fade_out_cycle)
    Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
    % Draw masks:
    Screen('DrawTextures', w, masksTex(fade_out_masks(frame)).all, [], this_trial__masks_position_xy, 180, [], fade_out_cycle(frame));
    Screen('Flip', w);
end


% End of block:
DrawFormattedText(w, 'Fin de bloque.', 'center', screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.right);
Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
Screen('Flip', w);


%% Save data
data.timestamps = timestamps;
data.trial_responses = trial_responses;
data.trial_response_times = trial_response_times;
data.frames_conds = frames_conds;
data.frames_faces = frames_faces;
data.frames_masks = frames_masks;
data.frames_alpha = 1-flipped;
data.baseline_framerate = baseline_framerate;
data.mask_framerate = mask_framerate;
data.oddball_framerate = oddball_framerate;



%% Close PsychToolbox

% Release keyboard with PsychToolbox, wait and clear workspace:
	KbQueueRelease();
	KbReleaseWait();
	WaitSecs(1);
	sca;
    clear Screen;
	ListenChar(0); % Restore keyboard output to Matlab


end


%% Set screen positions of fixed elements:
function [vergence_bars_positions, fix_position_left, fix_position_right, textBox, stereogram_box] =...
    create_static_images_boxes(expmnt, xCenter, yCenter, screenYpixels, mask_box)
    % create_static_images_boxes(expmnt, xCenter, yCenter, screenYpixels, mask_box)
    % 
    % Function to generate variables with screen positions of:
    % Fixation cross, vergence vars, and arrows.

    x_displacement_left =   expmnt.x_displacement * expmnt.cfs_1mirror_set_left_prop_of_displace;
    x_displacement_right =  expmnt.x_displacement * expmnt.cfs_1mirror_set_right_prop_of_displace;
    
    % =========================================================================
    % Fixation cross:
    % =========================================================================
    fix_position_left =  [xCenter - x_displacement_left,  yCenter]; 
    fix_position_right = [xCenter + x_displacement_right, yCenter]; 
    
    % =========================================================================
	% Vergence bar position:
    % =========================================================================
    if expmnt.cfs_with_one_mirror_set == 1
        vergence_bar_left =  expmnt.vergence_bar*expmnt.cfs_1mirror_set_size_ratio; % Original size: [0, 0, 16, 120]
        vergence_bar_right = expmnt.vergence_bar; % Original size: [0, 0, 16, 120]
        % This offplace is the same for both eyes, shouln't be, should be
        % smaller for left-eye when diff size: **********   PEND  **************
        vergence_bar_offplace_left =  expmnt.mask_width*expmnt.cfs_1mirror_set_size_ratio/2 + (20/(expmnt.mask_width*expmnt.cfs_1mirror_set_size_ratio))*expmnt.mask_width*expmnt.cfs_1mirror_set_size_ratio;
        vergence_bar_offplace_right = expmnt.mask_width/2 + (20/expmnt.mask_width)*expmnt.mask_width;
        vergence_bar_left =  CenterRectOnPointd(vergence_bar_left, xCenter, yCenter);
        vergence_bar_right = CenterRectOnPointd(vergence_bar_right, xCenter, yCenter);
        vergence_bar_left =     vergence_bar_left - [x_displacement_left, 0, x_displacement_left, 0];
        vergence_bar_right =    vergence_bar_right + [x_displacement_right, 0, x_displacement_right, 0];
        vergence_bars_positions = [ vergence_bar_left + [vergence_bar_offplace_left, 0, vergence_bar_offplace_left, 0];...
                                    vergence_bar_left - [vergence_bar_offplace_left, 0, vergence_bar_offplace_left, 0];...
                                    vergence_bar_right + [vergence_bar_offplace_right, 0, vergence_bar_offplace_right, 0];...
                                    vergence_bar_right - [vergence_bar_offplace_right, 0, vergence_bar_offplace_right, 0] ]';
    else
        vergence_bar_offplace = expmnt.mask_width/2 + 20;
		vergence_bar = CenterRectOnPointd(expmnt.vergence_bar, xCenter, yCenter);
		vergence_bar_right = vergence_bar + [expmnt.x_displacement, 0, expmnt.x_displacement, 0];
		vergence_bar_left = vergence_bar - [expmnt.x_displacement, 0, expmnt.x_displacement, 0];
		vergence_bars_positions = ...
			[ vergence_bar_right + [vergence_bar_offplace, 0, vergence_bar_offplace, 0];...
				vergence_bar_right - [vergence_bar_offplace, 0, vergence_bar_offplace, 0];...
				vergence_bar_left + [vergence_bar_offplace, 0, vergence_bar_offplace, 0];...
				vergence_bar_left - [vergence_bar_offplace, 0, vergence_bar_offplace, 0] ]';
    end
    
    % =========================================================================
    % Textboxes positions:
    % =========================================================================
    textBox_right = mask_box;
    textBox_left =  mask_box*expmnt.cfs_1mirror_set_size_ratio;
	textBox_right = CenterRectOnPointd(textBox_right, xCenter, yCenter);
    textBox_left =  CenterRectOnPointd(textBox_left, xCenter, yCenter);
	textBox.right = textBox_right + [x_displacement_right, 0, x_displacement_right, 0];
	textBox.left =  textBox_left - [x_displacement_left,  0, x_displacement_left, 0];
    
    
    % =========================================================================
    % Stereograms:
    % =========================================================================
    stereogram_box_left = CenterRectOnPointd(mask_box*expmnt.cfs_1mirror_set_size_ratio, xCenter, yCenter); 
    stereogram_box_right = CenterRectOnPointd(mask_box, xCenter, yCenter); 
    stereogram_box.left = stereogram_box_left - [x_displacement_left,  0, x_displacement_left, 0];
    stereogram_box.right = stereogram_box_right  + [x_displacement_right, 0, x_displacement_right, 0];
end


%% Load, edit and convert images into PsychToolbox textures:
function [imagesEdited, imagesReady, mask_frgd, mask_bkgd] = load_edit_convert_imgs(list_filenames_images, dir_images, expmnt, window)
    % load_edit_convert_imgs(list_filenames_images, dir_images, expmnt, window)
    %
    % LOAD:
    fprintf('\n\n============   LOADING IMAGES    ================\n');
    for i = 1:length(list_filenames_images)
        % Load images into Matlab separating fam and unfam faces:
        img_to_load     = list_filenames_images{i};
        fprintf('\n Load: %s  ', img_to_load);
        imagesLoaded{i} = imread(fullfile(dir_images, img_to_load));
        % Check if on grayscale, and convert to grayscale if so:
        if ndims(imagesLoaded{i}) == 3
            imagesLoaded{i} = rgb2gray(imagesLoaded{i});
        end
    end

    % EDIT: lumMatch, background colour, and convert into textures
    [y_size, x_size] = size(imagesLoaded{1});
    % Equalize images on sf, and then on hist, and set background to single color:
    mask_frgd = createEllipse(x_size, y_size);
    mask_bkgd = mask_frgd == 0;
    % Images should be already sfMatch & histMatch with 'exp06_edit_faces_SHINE':
    % SHINE receives imgs in uint8 (0-255) cell variables:
    imagesEdited =     lumMatch(imagesLoaded, mask_frgd, [expmnt.img_lum_mean, expmnt.img_contrast_sd]);
    % Set background colour to confirm it:
    for i = 1:length(imagesEdited)
        imagesEdited{i}(mask_bkgd) =   round(expmnt.background * expmnt.im_white); % Background colour
    end

    % If convolving with gaussian: to blur edges
    if expmnt.gaussian_envelope == 1
        x_radius_sd = 1.8; y_radius_sd = 1.8;
        x_size = size(imagesEdited{1}, 2);
        y_size = size(imagesEdited{1}, 1);
        msx = x_size/2; msy = y_size/2;
        [x, y] = meshgrid(-msx:msx-1, -msy:msy-1);
        xsd = msx/x_radius_sd; ysd = msy/y_radius_sd;
        mask_rect = uint8(255 - round(255 - exp(-((x/xsd).^2)-((y/ysd).^2))*255));
        for i = 1:length(imagesEdited)
            imagesEdited{i}(:,:,2) = mask_rect;
        end
    end

    % CONVERT INTO TEXTURES:
    for i = 1:length(imagesEdited)
        % Make the image into a texture:
        imagesReady(i).theTexture =      Screen('MakeTexture', window, imagesEdited{i});
    end
end
