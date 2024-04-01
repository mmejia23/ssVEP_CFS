function [data, expmnt] = ssvepCFS(subject, supps_alpha)
% Info
% 
% Pending:
% [X] Adjust mask and faces sizes: based on visual angle degrees.
% [X] Add selection of short trials: 10 famfaces cycles? Change in duration
% [ ] Add task: categorize Cneuro vs famous?
% [ ] Set parts mask alpha from UML estimation.

% For the setup:
% [ ] Code UML psychophysics estimation.
% [ ] What can work as *divider*?
% [ ] Do we need more unfamiliar faces?
% [ ] Do we need more/differente familiar faces?

clc;
fprintf(    '_________________________________________________________________________\n');
fprintf('\n\n      CFS & ssVEP experiment   \n\n');
fprintf(    '      Code: Manuel Mejia   \n');
fprintf(    '      Date: 2024-March-29   \n');
fprintf('\n\n_________________________________________________________________________\n');

addpath('SHINEtoolbox/');
addpath('convert_screen_sizes');
expmnt.savescreen =                 0;
expmnt.PTB_SkipSyncTexts_config =   1;
expmnt.debug_pt =                   0;
expmnt.BeampositionQueryWorkaround = 1;
expmnt.set_custom_screen_size =     0;

%% 60 or 85 Hz
% Frequencies that can be tagged depending on monitor screen rate:
% 1./([1:50].*(1000./60)).*1000 % 1.2, 6, 8.5714 = 50, 10, 7
% 1./([50, 10, 7].*(1000./60)).*1000
% If we define Hz for masks and baseline stimuli:
% 1./([1:80].*(1000./85)).*1000 % 1.2143, 6.0714, 8.5  = 70, 14, 10
% 1./([70, 14, 10].*(1000./85)).*1000

%% Set experiment variables
%__________________________________________________________________________
dir_images = '../stimuli/caras_experimento/';
filename_structure = '../data/sub-%s_task-ssvepCFS';
expmnt.dir_of_texts_png_files = '../stimuli/textos/';

% Blocks and trials
expmnt.trial_duration =     10;
expmnt.block_conds_order =  [0, 1, 0, 1];

% Dynamics
expmnt.masks_hz =       8.5; % approx, is defined according to monitor framerate
expmnt.baseline_hz =    6;
expmnt.fade_in_time =   2;
expmnt.fade_out_time =  2;
expmnt.odd_frequency =  5; % each 5th baseline stim, the oddball

% Suppressions: based on averages of pupillometry study
% expmnt.default_alpha_masks_supp1   = 0.04;
% expmnt.default_alpha_masks_supp2   = 0.32;
% expmnt.default_alpha_masks_supp3   = 0.60;
if nargin > 1
    expmnt.supp1__mask_alpha =  supps_alpha(1);
    expmnt.supp2__mask_alpha =  supps_alpha(2);
else
    expmnt.supp1__mask_alpha =  0.03;
    expmnt.supp2__mask_alpha =  0.32;
end

% Sizes

% Stimuli:
expmnt.stimuli_width =          100;
expmnt.stimuli_height =         133;

% Masks box (same as stereograms and text boxes):
expmnt.mask_width =             240;
expmnt.mask_height =            240;
expmnt.num_masks =              500; %approx

% Vergence bars
expmnt.vergence_bar_width =     32;
expmnt.vergence_bar_height =    240;

expmnt.im_white    = 255;
expmnt.img_contrast_sd         = 12; % approx 20% Michelson contrast. Need to check. %%%  EDIT BEFORE EXP. %%% 
expmnt.cfs_masks_contrast_sd   = 82; % the mean of the original function make_mondrian_masks
expmnt.gaussian_envelope   = 0; % If blur edges with gaussian envelope.

expmnt.position_task_offset    = 0;
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


%% Filename
%__________________________________________________________________________
filename = sprintf(filename_structure, subject);
if exist([filename, '_beh.mat'], 'file')
    % If there are previous files, rename them to separate them:
    matFilesFound           = dir([filename '*']);
    matFilesFound           = {matFilesFound.name};
    newname_prev_matFile    = [filename, '_desc-try', num2str(length(matFilesFound)), '_beh.mat'];
    movefile([filename, '_beh.mat'], newname_prev_matFile);
end


%% Estimated vars within experiment:
%__________________________________________________________________________

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
%__________________________________________________________________________

% Get list of filenames of stimuli:
familiar_images = dir([dir_images filesep 'familiar*.png']);
familiar_images = {familiar_images.name};
unfamiliar_images = dir([dir_images filesep 'unfamiliar*.png']);
unfamiliar_images = {unfamiliar_images.name};
assert(length(familiar_images)>0 & length(unfamiliar_images)>0);


%% Start PsychToolbox
%__________________________________________________________________________

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
%__________________________________________________________________________
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
%__________________________________________________________________________
% Open an on screen window:
[w, wRect] = PsychImaging('OpenWindow', screenNumber, expmnt.background);
% Get size of on screen window:
[screenXpixels, screenYpixels] = Screen('WindowSize', w);
% Query the frame duration. It will be used to control the while loop.
ifi = Screen('GetFlipInterval', w);
hz = Screen('FrameRate', w); % Nominal frame rate
% Get the centre coordinate of the window.
[xCenter, yCenter] = RectCenter(wRect);
% Set up alpha-blending for smooth (anti-aliased) lines????
Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('HideCursorHelper', w);


%% Estimate stimuli, masks, vergence bars, fixation, and texts sizes in px
%__________________________________________________________________________
screen_size_px =    [screenXpixels, screenYpixels]; 	% width: 800px, height: 600px
screen_distance_mm = [635]; 	% eye to center of screen distance: 635 mm.
screen_size_mm =    [330, 220]; 	% width: ?mm, height: ?mm.
stim_center_deg =   [0, 0]; 		% The center of the stimuli matches the center of the screen: degrees of eccentricity

% Stimuli:
stim_size_deg =     [3.1, 4.1]; 		% width: , height: 
[stim_size_px, stim_size_mm, stim_center_px] =...
		deg2px2(stim_size_deg, screen_size_px, screen_distance_mm, screen_size_mm, stim_center_deg);

% Masks/stereograms/text boxes:
mask_size_deg =     [4.5, 4.5]; 		% width: , height: 
[mask_size_px, mask_size_mm, mask_center_px] =...
		deg2px2(mask_size_deg, screen_size_px, screen_distance_mm, screen_size_mm, stim_center_deg);

% Vergence bars:
vergence_size_deg =     [0.6, 4.5]; 		% width: , height: 
[vergence_size_px, vergence_size_mm, vergence_center_px] =...
		deg2px2(vergence_size_deg, screen_size_px, screen_distance_mm, screen_size_mm, stim_center_deg);

% Set sizes:
expmnt.vergence_bar_width =     round(vergence_size_px(1));
expmnt.vergence_bar_height =    round(vergence_size_px(2));
expmnt.mask_width =             round(mask_size_px(1));
expmnt.mask_height =             round(mask_size_px(2));
expmnt.stimuli_width =          round(stim_size_px(1));
expmnt.stimuli_height =         round(stim_size_px(2));

% Target box:
target_box = [0 0 expmnt.stimuli_width expmnt.stimuli_height]; % Position of target.
% Masks box:
mask_box = [0, 0, expmnt.mask_width, expmnt.mask_height];
expmnt.vergence_bar = [0, 0, expmnt.vergence_bar_width, expmnt.vergence_bar_height];

%% Prepare conditions vectors
%__________________________________________________________________________

expmnt.monitor_hz =     hz;
mask_framerate =        round(expmnt.monitor_hz/expmnt.masks_hz);
baseline_framerate =    round(expmnt.monitor_hz/expmnt.baseline_hz);
oddball_framerate =     baseline_framerate*expmnt.odd_frequency;

% Just a vector with ones for each frame for an entire block:
frames_conds = ones(1, expmnt.trial_duration*expmnt.monitor_hz);


%% Duty cycle for sinusoidal modulation
%__________________________________________________________________________
flickerRate =   expmnt.baseline_hz; % 6 Hz
cycle =         round(linspace(0, 360, hz));  %round(1/(frames)*1000)));
cycle =         repmat(cycle,[1,expmnt.trial_duration,1]);
trans =         cosd(cycle*flickerRate+180)*0.5+0.5; %Create sinusoidal levels for pixel level
% trans(end) = [];                                    %Remove last element (as first is the same)
%180 = shift phase; *0.5+0.5 = Transparency level 0-1
%mins = findminima(trans);
flipped = 1 - trans;
[peaks,locs]=findpeaks(flipped);      %Flip then find troughs
mins = locs;


%% ------------------------------------------------------------------------
% Load vergence bar:
%__________________________________________________________________________
vergence(1).img = imread(expmnt.vergence_bar_file);
vergence(1).img = mat2gray(vergence(1).img, [0, 3]);
if expmnt.show_vergence_bars == 0
    vergence(1).img(:, :, 2) = zeros(size(vergence(1).img, 1), size(vergence(1).img, 2));
end
vergence(1).tex = Screen('MakeTexture', w, vergence(1).img);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load stereograms:
%__________________________________________________________________________
stereograms_list = dir(fullfile(expmnt.dir_of_stereograms, '*l.jpg'));
for i = 1:length(stereograms_list)
    stereograms(i).img_left =  imread(fullfile(expmnt.dir_of_stereograms, sprintf('stereogram%02d_l.jpg', i)));
    stereograms(i).img_right = imread(fullfile(expmnt.dir_of_stereograms, sprintf('stereogram%02d_r.jpg', i)));

    stereograms(i).tex_left =  Screen('MakeTexture', w, stereograms(i).img_left);
    stereograms(i).tex_right = Screen('MakeTexture', w, stereograms(i).img_right);
end
clearvars stereograms_list;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load text images:
%__________________________________________________________________________
% Get list of filenames of texts:
texts_directory = dir(fullfile(expmnt.dir_of_texts_png_files, 'texto*'));
% Get list of images:
texts_images = {texts_directory.name};
for i = 1:length(texts_images)
    [temp, ~, alpha] = imread(fullfile(expmnt.dir_of_texts_png_files, texts_images{i}));
    textsImages(i).images = cat(3, temp, temp, temp, alpha);
    textsImages(i).textures = Screen('MakeTexture', w, textsImages(i).images);
end
clearvars texts_directory;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimuli fixed positions:
%__________________________________________________________________________
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
%__________________________________________________________________________

Screen('TextSize', w, font_size);
KbQueueCreate([], responseKeys);
KbQueueStart;
KbQueueFlush();
x_displacement_set = 0;
while x_displacement_set == 0
    % Draw stereograms:
    Screen('DrawTextures', w, stereograms(which_stereogram(1)).tex_left, [], stereogram_box.left, []);
    Screen('DrawTextures', w, stereograms(which_stereogram(1)).tex_right, [], stereogram_box.right, []);
    Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
    Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
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
WaitSecs(.1);

% Preparando imagenes...
Screen('DrawTextures', w, textsImages(1).textures, [], [textBox.left; textBox.right]', 0, [], []);
%     DrawFormattedText(w, 'Preparando las imagenes...', 'center', screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.left);
%     DrawFormattedText(w, 'Preparando las imagenes...', 'center', screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.right);
Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
Screen('Flip', w);
WaitSecs(.1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and edit stimuli for this participant:
%__________________________________________________________________________

% Edit images: lumMatch, background colour, and convert into textures
% All faces
if 1
    [loadedFacesEdited, loadedFaces, mask_frgd, mask_bkgd] =...
        load_edit_convert_imgs([unfamiliar_images, familiar_images], dir_images, expmnt, w);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create and prepare masks:
%__________________________________________________________________________
% make_mondrian_masks(sz_x, sz_y, n_masks, shape, selection)
% shape=2, circles. % selection=2, grayscale.
expmnt.num_masks =  ceil(length(frames_conds)/mask_framerate); % Minimum
temp =              Expand(1:expmnt.num_masks, mask_framerate, 1);
frames_masks =      temp(1:length(frames_conds));
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
%__________________________________________________________________________
Screen('DrawTextures', w, textsImages(2).textures, [], [textBox.left; textBox.right]', 0, [], []);
% DrawFormattedText(w, 'Imagenes listas.', 'center', screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.right);
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
%__________________________________________________________________________
% Draw stereograms:
Screen('DrawTextures', w, stereograms(which_stereogram(2)).tex_left, [], stereogram_box.left, []);
Screen('DrawTextures', w, stereograms(which_stereogram(2)).tex_right, [], stereogram_box.right, []);
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


for block = 1:length(expmnt.block_conds_order)

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % screen: press any key to begin. 'Presiona cualquier tecla \n\npara iniciar.'
    %__________________________________________________________________________
    Screen('DrawTextures', w, textsImages(3).textures, [], [textBox.left; textBox.right]', 0, [], []);
    % DrawFormattedText(w, 'Presione cualquier tecla para iniciar.', 'center',...
    %     screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.right);
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
    
    %% Start block loop vars:
    %__________________________________________________________________________
    timestamps =            nan(1, length(frames_conds));
    trial_responses =       [];
    trial_response_times =  [];
    trial_response =        [];
    trial_response_time =   [];
    
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
    temp = Shuffle((1:length(familiar_images))+length(unfamiliar_images));
    frames_faces(idx_oddball) = temp(1:length(frames_faces(idx_oddball)));
    clearvars temp;

    % Expand each image index to the number of frames to be presented:
    frames_faces = Expand(frames_faces, baseline_framerate, 1);
    
    % Define mask alpha: 0=visible, 1=suppressed
    switch expmnt.block_conds_order(block)
        case 0
            this_trial__mask_alpha = expmnt.supp1__mask_alpha;
        case 1
            this_trial__mask_alpha = expmnt.supp2__mask_alpha;
    end
    
    % Mask position and potency:
    this_trial__xy_target_masks_displacement = ...
                [expmnt.x_displacement, 0, expmnt.x_displacement, 0];
    this_trial__masks_position_xy =...
        mask_box_right + expmnt.cfs_1mirror_set_right_prop_of_displace * this_trial__xy_target_masks_displacement;

    % Target position:
    this_trial__target_position_xy = ...
        target_box_left - expmnt.cfs_1mirror_set_right_prop_of_displace * this_trial__xy_target_masks_displacement;

    %% Block loop
    %__________________________________________________________________________
    KbQueueFlush();

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
        Screen('DrawTextures', w, masksTex(fade_in_masks(frame)).all, [],...
            this_trial__masks_position_xy, 180, [], fade_in_cycle(frame));
        Screen('Flip', w);
    end

    for frame = 1:length(frames_conds)

        Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
        % Draw target:
        Screen('DrawTextures', w, loadedFaces(frames_faces(frame)).theTexture, [],...
            this_trial__target_position_xy, 0, [], 1-flipped(frame));
        % Draw masks:
        Screen('DrawTextures', w, masksTex(frames_masks(frame)).all, [],...
            this_trial__masks_position_xy, 180, [], this_trial__mask_alpha);

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


    %% End of block:
    %__________________________________________________________________________
    Screen('DrawTextures', w, textsImages(6).textures, [], [textBox.left; textBox.right]', 0, [], []);
    % DrawFormattedText(w, 'Fin de bloque.', 'center', screenYpixels * pos_first_text_line, [1 1 1], [], [], [], [], [], textBox.right);
    Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
    Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
    Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
    Screen('Flip', w);
    KbWait;
    WaitSecs(.5);


    %% Save data
    %__________________________________________________________________________
    data(block).subject =               subject;
    data(block).timestamps =            timestamps;
    data(block).trial_responses =       trial_responses;
    data(block).trial_response_times =  trial_response_times;
    data(block).frames_conds =          frames_conds;
    data(block).frames_faces =          frames_faces;
    data(block).frames_masks =          frames_masks;
    
    expmnt.frames_alpha =           1-flipped;
    expmnt.baseline_framerate =     baseline_framerate;
    expmnt.mask_framerate =         mask_framerate;
    expmnt.oddball_framerate =      oddball_framerate;
    expmnt.onscreen_positions.vergence_bars_positions =     vergence_bars_positions;
    expmnt.onscreen_positions.fix_position_left =           fix_position_left;
    expmnt.onscreen_positions.fix_position_right =          fix_position_right;
    expmnt.onscreen_positions.textBox =                     textBox;
    expmnt.onscreen_positions.stereogram_box =              stereogram_box;
    expmnt.onscreen_positions.masks_position_xy =           this_trial__masks_position_xy;
    expmnt.onscreen_positions.target_position_xy =          this_trial__target_position_xy;
    expmnt.onscreen_positions.ifi =                         ifi;
    expmnt.onscreen_positions.screenXpixels =               screenXpixels;
    expmnt.onscreen_positions.screenYpixels =               screenYpixels;

    save([filename, '_beh.mat'], 'data', 'expmnt');
    
end

%% Close PsychToolbox
%__________________________________________________________________________

% Release keyboard with PsychToolbox, wait and clear workspace:
KbQueueRelease();
KbReleaseWait();
WaitSecs(1);
sca;
clear Screen;
ListenChar(0); % Restore keyboard output to Matlab

%% Plot test of cycles:
%__________________________________________________________________________
figure;
subplot(4,1,1);
plot(data(1).frames_faces<73); 
title(sprintf('Unfamiliar faces: %.02f Hz', expmnt.baseline_hz));
ylim([-0.1, 1.1]);
subplot(4,1,2);
plot(data(1).frames_faces>72); 
title(sprintf('Familiar faces: %.02f Hz', expmnt.baseline_hz/expmnt.odd_frequency));
ylim([-0.1, 1.1]);
subplot(4,1,3);
plot(diff(data(1).frames_masks));
title(sprintf('CFS Masks: %.02f Hz', expmnt.masks_hz));
ylim([-0.1, 1.1]);
subplot(4,1,4);
plot(diff(data(1).timestamps));
title(sprintf('IFI, should be around: %.03f s', expmnt.onscreen_positions.ifi));

% figure;
% plot(data.frames_faces);
% title('Unfamiliar faces: 1-72, Familiar faces: 73-144');

end


%% Set screen positions of fixed elements:
%__________________________________________________________________________
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
        vergence_bar_offplace_left =  expmnt.mask_width*expmnt.cfs_1mirror_set_size_ratio/2 +...
            ((expmnt.vergence_bar_width*1.25/2)/(expmnt.mask_width*expmnt.cfs_1mirror_set_size_ratio))*...
            expmnt.mask_width*expmnt.cfs_1mirror_set_size_ratio;
        vergence_bar_offplace_right = expmnt.mask_width/2 + ...
            ((expmnt.vergence_bar_width*1.25/2)/expmnt.mask_width)*expmnt.mask_width;
        vergence_bar_left =  CenterRectOnPointd(vergence_bar_left, xCenter, yCenter);
        vergence_bar_right = CenterRectOnPointd(vergence_bar_right, xCenter, yCenter);
        vergence_bar_left =     vergence_bar_left - [x_displacement_left, 0, x_displacement_left, 0];
        vergence_bar_right =    vergence_bar_right + [x_displacement_right, 0, x_displacement_right, 0];
        vergence_bars_positions = [ vergence_bar_left + [vergence_bar_offplace_left, 0, vergence_bar_offplace_left, 0];...
                                    vergence_bar_left - [vergence_bar_offplace_left, 0, vergence_bar_offplace_left, 0];...
                                    vergence_bar_right + [vergence_bar_offplace_right, 0, vergence_bar_offplace_right, 0];...
                                    vergence_bar_right - [vergence_bar_offplace_right, 0, vergence_bar_offplace_right, 0] ]';
    else
        vergence_bar_offplace = expmnt.mask_width/2 + (expmnt.vergence_bar_width*1.25/2);
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
%__________________________________________________________________________
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

%% Done:
%__________________________________________________________________________
% [X] Custom adjustment of x_displacement.
% [X] Present masks depending on which frame we are.
% [X] Present stimuli depending on which frame we are.
% [X] Add fade-in & fade-out for each stimulus.
% [X] Test keyboard check (do we have problems with timing?).
%     Apparently there is no lag!!
% [X] Save to Github.
% [X] Add overall fade-in & fade-out for entire block.
% [X] Add values for mask alpha levels.
% [X] Verify frames presented: save conditions vectors.
% [X] Add stereogram at first screen for position calibration.
% [X] Save mask, target, vergence bars positions: expmnt.
% [X] Change place of Hz estimations (after PsychToolbox start).
% [X] Save data from each block at the end of each block. BIDS format.
% [X] Add selection of visible and invisible blocks.
