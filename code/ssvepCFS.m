function [data, expmnt] = ssvepCFS(subject, varargin)
% Info
% 
% Pending:
% [X] Add four conditions: visible vs invisible, upright vs inverted faces
% [ ] Prepare at least 15 photos (natural) of 10 identities.
% [ ] Show screen with face exemplars for response.
% [X] Clean response collection, missing some trials.. 1 and 2 each block
% [ ] Confirm that PF task works well and saves data correctly.
% [ ] Prepare way to use UML object to define part's thresholds.
% [ ] Estimate times: PF 5 mins approx? Main task 15 mins approx?
% [ ] Sync with triggerbox.
% [X] Add more stereograms and show different one each block.

% [ ] Do we need more unfamiliar faces? Probably yes: use webmorphR
% [ ] Do we need more/different familiar faces? Probably yes: use webmorphR
% [ ] Fotos completas
% [X] Que la ultima no sea familiar
% [ ] Agregar la cantidad de estimulos para igualar a Rossion
% [ ] Programar condicion con un solo estimulo central
% [ ] Programar condicion sin CFS


clc;
fprintf(    '_________________________________________________________________________\n');
fprintf('\n\n      CFS & ssVEP experiment   \n\n');
fprintf(    '      Code: Manuel Mejia   \n');
fprintf(    '      Date: 2024-March-29   \n');
fprintf('\n\n_________________________________________________________________________\n');

addpath('SHINEtoolbox/');
addpath('UML_ver2.2/');
addpath('convert_screen_sizes');
expmnt.savescreen =                 0;
expmnt.PTB_SkipSyncTexts_config =   1;
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
dir_images_rep = 'C:\Users\cthul\Dropbox\PROYECTOS\0-PhD\ssVEP_CFS\stimuli\caras_familiares_edited';
filename_structure = '../data/sub-%s_task-%s';
expmnt.dir_of_texts_png_files = '../stimuli/textos/';

% Dynamics
expmnt.masks_hz =       8.5; % approx, is defined according to monitor framerate
expmnt.baseline_hz =    3; % We can tinker this one; 3; 6;
expmnt.fade_in_time =   2;
expmnt.fade_out_time =  2;
expmnt.odd_frequency =  5; % each 5th baseline stim, the oddball

% Suppressions: based on averages of pupillometry study
% expmnt.default_alpha_masks_supp1   = 0.04;
% expmnt.default_alpha_masks_supp2   = 0.32;
% expmnt.default_alpha_masks_supp3   = 0.60;
if any(find(strcmp(varargin, 'alphas')))
    expmnt.supp1__mask_alpha =  varargin{find(strcmp(varargin, 'alphas')) + 1}(1);
    expmnt.supp2__mask_alpha =  varargin{find(strcmp(varargin, 'alphas')) + 1}(2);
else
    expmnt.supp1__mask_alpha =  0.03;
    expmnt.supp2__mask_alpha =  0.32;
end

if any(find(strcmp(varargin, 'pf')))
    pf_estimation = varargin{find(strcmp(varargin, 'pf')) + 1};
else
    pf_estimation = 0;
end

if any(find(strcmp(varargin, 'debug')))
    expmnt.debug_pt = varargin{find(strcmp(varargin, 'debug')) + 1};
else
    expmnt.debug_pt = 0;
end

% Blocks and trials
if pf_estimation
    expmnt.stim_per_trial =     1;
    expmnt.trial_duration =     1/expmnt.baseline_hz;
    expmnt.n_trials_per_block = 5;
    expmnt.block_conds_order =  [0, 0, 0, 0];
    % PF settings
    expmnt.constant_stimuli_method =    0;
    expmnt.pf_type_estimations =        1;
    expmnt.pf_trials_each_per_block =   35;
else
    expmnt.stim_per_trial =     expmnt.odd_frequency * 2 + 3;
    expmnt.trial_duration =     expmnt.stim_per_trial ./ expmnt.baseline_hz;
    expmnt.n_trials_per_block = 3;
    expmnt.block_conds_order =  [2, 3];
end
expmnt.trial_conds_order =  repmat([0,1], length(expmnt.block_conds_order), expmnt.n_trials_per_block);
expmnt.identity_orders = {'BarackObama','RaulCastro','DiazCanel','Beyonce','OmaraPortuondo'};
expmnt.identity_orders = repmat(expmnt.identity_orders,...
            length(expmnt.block_conds_order),...
            ceil(expmnt.n_trials_per_block./length(expmnt.identity_orders)));
        
%% Sizes

% Screen sizes:
% Eyelink lab: distance = 635mm, screen size = 330x220 mm.
screen_distance_mm =    [635]; 	% eye to center of screen distance: ?? mm.
screen_size_mm =        [520, 290]; 	% width: ?mm, height: ?mm.
% EEG lab with video beam: distance = 700mm?, screen size = 1520x950mm
% screen_distance_mm =    [1000]; 	% eye to center of screen distance: ?? mm.
% screen_size_mm =        [1520, 950]; 	% width: ?mm, height: ?mm.
expmnt.custom_screen_center =   0; % 1 for EEG clinic at CNeuro.
expmnt.xCenter_custom =         795;
expmnt.yCenter_custom =         520;

% Stimuli:
% stim_size_deg =     [3.1, 4.1]; 		% width: , height:
stim_size_deg =     [5, 6.6130]; 
% Set this default, but later is changed:
expmnt.stimuli_width =          100;
expmnt.stimuli_height =         133;

% Masks box (same as stereograms and text boxes):
% mask_size_deg =     [4.5, 4.5]; 		% width: , height: 
mask_size_deg =     [5.5, 5.5];
% Set this default, but later is changed:
expmnt.mask_width =             240;
expmnt.mask_height =            240;
expmnt.num_masks =              500; %approx

% Vergence bars
% vergence_size_deg =     [0.6, 4.5]; 		% width: , height: 
vergence_size_deg =     [0.73, 5.5];
% Set this default, but later is changed:
expmnt.vergence_bar_width =     32;
expmnt.vergence_bar_height =    240;

expmnt.im_white    = 255;
expmnt.img_contrast_sd         = 12; % approx 20% Michelson contrast. Need to check. %%%  EDIT BEFORE EXP. %%% 
expmnt.cfs_masks_contrast_sd   = 82; % the mean of the original function make_mondrian_masks
expmnt.gaussian_envelope   = 0; % If blur edges with gaussian envelope.

expmnt.position_task_offset =   0;
expmnt.show_vergence_bars =     1;
expmnt.vergence_bar_file =      '../stimuli/texture_vertical.bmp';
expmnt.dir_of_stereograms =     '../stimuli/stereograms2';
expmnt.stereogram_shapes =      {'square', 'T', 'square', 'square', 'upper-square', 'lower-square', 'two-squares', 'L', 'U', 'big-square'};
expmnt.stereogram_for_pf               = [1, 2, 3, 8, 5, 9, 10, 13, 15];
expmnt.stereogram_for_exp              = [1, 2, 3, 8, 5, 9, 10, 13, 15];

expmnt.fix_color = [1,1,1];
expmnt.img_lum_mean    = 120; % background x im_white, equal to background mean luminance.
expmnt.background = 0.47;
expmnt.cfs_with_one_mirror_set = 0; % Set to 1 if only using mirrors for one eye: for using Eyelink tracking.
expmnt.x_displacement = 300;
expmnt.y_displacement = 0;
if expmnt.cfs_with_one_mirror_set == 1
    expmnt.cfs_1mirror_set_left_prop_of_displace       = 0.90;
    expmnt.cfs_1mirror_set_right_prop_of_displace      = 1.10;
    expmnt.cfs_1mirror_set_size_ratio                  = 0.86;
else
    expmnt.cfs_1mirror_set_left_prop_of_displace       = 1.00;
    expmnt.cfs_1mirror_set_right_prop_of_displace      = 1.00;
    expmnt.cfs_1mirror_set_size_ratio                  = 1.00;
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
block_font_size = 15;
    
%% Prepare stimuli
%__________________________________________________________________________

% Get list of filenames of stimuli:
familiar_faces = dir([dir_images filesep 'familiar*.png']);
familiar_faces = {familiar_faces.name};
unfamiliar_faces = dir([dir_images filesep 'unfamiliar*.png']);
unfamiliar_faces = {unfamiliar_faces.name};
unfamiliar_houses = dir([dir_images filesep 'unfamiliar_house*.png']);
unfamiliar_houses = {unfamiliar_houses.name};
familiar_faces_rep = dir([dir_images_rep filesep 'familiar*.png']);
familiar_faces_rep = {familiar_faces_rep.name};
assert(length(familiar_faces)>0 & ...
    length(unfamiliar_faces)>0 & ...
    length(unfamiliar_houses)>0 & ...
    length(familiar_faces_rep)>0);


%% Filename
%__________________________________________________________________________
if pf_estimation == 1
    filename = sprintf(filename_structure, subject, 'pf');
else
    filename = sprintf(filename_structure, subject, 'ssvepCFS');
end
if exist([filename, '_beh.mat'], 'file')
    % If there are previous files, rename them to separate them:
    matFilesFound           = dir([filename '*']);
    matFilesFound           = {matFilesFound.name};
    newname_prev_matFile    = [filename, '_desc-try', num2str(length(matFilesFound)), '_beh.mat'];
    movefile([filename, '_beh.mat'], newname_prev_matFile);
end


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
fix_xCoords =   [-fix_size, fix_size, 0,           0];
fix_yCoords =   [0,         0,        -fix_size,   fix_size];
fixCoords =         [fix_xCoords; fix_yCoords];
fixCoordsLeft =     fixCoords*expmnt.cfs_1mirror_set_size_ratio;
fixCoordsRight =    fixCoords;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set available keys to press:
%__________________________________________________________________________
% Get key codes:
escapeKey =     KbName('ESCAPE');
spaceKey =      KbName('space');
keyUp =          KbName('UpArrow');
keyDown =          KbName('DownArrow');
keyLeft =          KbName('LeftArrow');
keyRight =          KbName('RightArrow');
key1 =          KbName('1!');
key2 =          KbName('2@');
tailorKeys =  zeros(1, 256);
tailorKeys([escapeKey, spaceKey, keyUp, keyDown, keyLeft, keyRight, key1, key2]) = 1;
responseKeys =  zeros(1, 256);
responseKeys([escapeKey, spaceKey, keyUp, keyDown, keyLeft, keyRight]) = 1;
pauseKeys =     zeros(1, 256);
pauseKeys(spaceKey) = 1;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start PsychToolbox window:
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
if expmnt.custom_screen_center
    xCenter = expmnt.xCenter_custom;
    yCenter = expmnt.yCenter_custom;
end
% Set up alpha-blending for smooth (anti-aliased) lines????
Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% Screen('HideCursorHelper', w);


%% Estimate stimuli, masks, vergence bars, fixation, and texts sizes in px
%__________________________________________________________________________
screen_size_px =    [screenXpixels, screenYpixels]; 	% width: 800px, height: 600px
stim_center_deg =   [0, 0]; 		% The center of the stimuli matches the center of the screen: degrees of eccentricity

% Stimuli:
% stim_size_deg =     [3.1, 4.1]; 		% width: , height:
stim_size_deg =     [5, 6.6130]; 
[stim_size_px, stim_size_mm, stim_center_px] =...
		deg2px2(stim_size_deg, screen_size_px, screen_distance_mm, screen_size_mm, stim_center_deg);

% Masks/stereograms/text boxes:
% mask_size_deg =     [4.5, 4.5]; 		% width: , height: 
mask_size_deg =     [5.5, 5.5];
[mask_size_px, mask_size_mm, mask_center_px] =...
		deg2px2(mask_size_deg, screen_size_px, screen_distance_mm, screen_size_mm, stim_center_deg);

% Vergence bars:
% vergence_size_deg =     [0.6, 4.5]; 		% width: , height: 
vergence_size_deg =     [0.73, 5.5];
[vergence_size_px, vergence_size_mm, vergence_center_px] =...
		deg2px2(vergence_size_deg, screen_size_px, screen_distance_mm, screen_size_mm, stim_center_deg);
    
% Response screen items:
responseItem_size_deg =     mask_size_deg./4;
[responseItem_size_px, responseItem_size_mm, responseItem_center_px] =...
		deg2px2(responseItem_size_deg, screen_size_px, screen_distance_mm, screen_size_mm, stim_center_deg);

% Set sizes:
expmnt.vergence_bar_width =     round(vergence_size_px(1));
expmnt.vergence_bar_height =    round(vergence_size_px(2));
expmnt.mask_width =             round(mask_size_px(1));
expmnt.mask_height =            round(mask_size_px(2));
expmnt.stimuli_width =          round(stim_size_px(1));
expmnt.stimuli_height =         round(stim_size_px(2));
expmnt.x_displacement =         round(expmnt.mask_width/2 + expmnt.vergence_bar_width*4);

% Target box:
target_box =    [0 0 expmnt.stimuli_width expmnt.stimuli_height]; % Position of target.
% Masks box:
mask_box =      [0, 0, expmnt.mask_width, expmnt.mask_height];
expmnt.vergence_bar = [0, 0, expmnt.vergence_bar_width, expmnt.vergence_bar_height];
responseItem_box = [0, 0, responseItem_size_px(1), responseItem_size_px(2)];

%% Prepare conditions vectors
%__________________________________________________________________________

expmnt.monitor_hz =     hz;
mask_framerate =        round(expmnt.monitor_hz/expmnt.masks_hz);
baseline_framerate =    round(expmnt.monitor_hz/expmnt.baseline_hz);
oddball_framerate =     baseline_framerate*expmnt.odd_frequency;

expmnt.masks_hz =       expmnt.monitor_hz / mask_framerate; % Real Hz
expmnt.baseline_hz =    expmnt.monitor_hz / baseline_framerate;

n_frames_per_trial = expmnt.stim_per_trial * baseline_framerate;

% Just a vector with ones for each frame for an entire block:
frames_conds = ones(1, n_frames_per_trial);


%% Duty cycle for sinusoidal modulation
%__________________________________________________________________________
flickerRate =   expmnt.baseline_hz; % 6 Hz
cycle =         round(linspace(0, 360, hz));  %round(1/(frames)*1000)));
cycle =         repmat(cycle,[1,ceil(expmnt.trial_duration),1]);
trans =         cosd(cycle*flickerRate+180)*0.5+0.5; %Create sinusoidal levels for pixel level
% trans(end) = [];                                    %Remove last element (as first is the same)
%180 = shift phase; *0.5+0.5 = Transparency level 0-1
%mins = findminima(trans);
flipped = 1 - trans;
[peaks,locs]=findpeaks(flipped);      %Flip then find troughs
mins = locs;
% Cut flipped vector for the length of the trial:
flipped = flipped(1:n_frames_per_trial);

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
left_stereograms_list = dir(fullfile(expmnt.dir_of_stereograms, '*l.png'));
right_stereograms_list = dir(fullfile(expmnt.dir_of_stereograms, '*r.png'));
for i = 1:length(left_stereograms_list)
    stereograms(i).img_left =  imread(fullfile(expmnt.dir_of_stereograms, left_stereograms_list(i).name));
    stereograms(i).img_right = imread(fullfile(expmnt.dir_of_stereograms, right_stereograms_list(i).name));

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



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First screen: loading images. Text: 'Preparando las imagenes...'
%__________________________________________________________________________

Screen('TextSize', w, font_size);
KbQueueCreate([], tailorKeys);
KbQueueStart;
KbQueueFlush();
x_displacement_set = 0;
while x_displacement_set == 0
    [vergence_bars_positions, fix_position_left, fix_position_right, textBox, stereogram_box] =...
            create_static_images_boxes(expmnt, xCenter, yCenter, screenYpixels, mask_box);
    % Draw stereograms:
    Screen('DrawTextures', w, stereograms(which_stereogram(1)).tex_left, [], stereogram_box.left, []);
    Screen('DrawTextures', w, stereograms(which_stereogram(1)).tex_right, [], stereogram_box.right, []);
%     Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
%     Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
    Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
    Screen('Flip', w);
    [pressed, firstPress] = KbQueueCheck;
    if pressed
        if firstPress(keyLeft)
            expmnt.x_displacement = expmnt.x_displacement - 5;
        end
        if firstPress(keyRight)
            expmnt.x_displacement = expmnt.x_displacement + 5;
        end
        if firstPress(keyUp)
            yCenter = yCenter - 5;
        end
        if firstPress(keyDown)
            yCenter = yCenter + 5;
        end
        if firstPress(spaceKey)
            x_displacement_set = 1;
        end
        if firstPress(key1)
            xCenter = xCenter - 5;
        end
        if firstPress(key2)
            xCenter = xCenter + 5;
        end
        if firstPress(escapeKey)
            closePT
            return
        end
    end
end
expmnt.xCenter = xCenter;
expmnt.yCenter = yCenter;
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stimuli fixed positions:
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

% Response screen positions:
responseItem_box = CenterRectOnPointd(responseItem_box, xCenter, yCenter);
responseItem_box = responseItem_box + ...
    [expmnt.x_displacement, expmnt.y_displacement,...
    expmnt.x_displacement, expmnt.y_displacement];
% 9 response options, in a 3x3 grid:
responseItem_rel_pos = [-.25,.25,-.25,.25; 0,.25,0,.25; .25,.25,.25,.25; ...
                        -.25,0,-.25,0; 0,0,0,0; .25,0,.25,0;...
                        -.25,-.25,-.25,-.25; 0,-.25,0,-.25; .25,-.25,.25,-.25];
responseItem_rel_pos = responseItem_rel_pos .* expmnt.mask_width;
responseItem_pos_all = repmat(responseItem_box, 9, 1);
responseItem_pos_all = responseItem_pos_all + responseItem_rel_pos;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and edit stimuli for this participant:
%__________________________________________________________________________

% Edit images: lumMatch, background colour, and convert into textures
% All faces
if 1
    [loadedStimEdited, loadedStim, mask_frgd, mask_bkgd] =...
        load_edit_convert_imgs([unfamiliar_faces, familiar_faces, unfamiliar_houses],...
        dir_images, expmnt, w);
end
if 1
    [loadedFamEdited, loadedFam, ~, ~] =...
        load_edit_convert_imgs([familiar_faces_rep],...
        dir_images_rep, expmnt, w);
end
% Define structure with characteristics of each image:
images_cell = [unfamiliar_faces, familiar_faces, unfamiliar_houses];
images_str = regexp(images_cell,...
    '(?<cond>\w+)_(?<stim>\w+)_(?<where>\w+)_(?<sex>\w+)_(?<name>\w+)_(?<num>\w+).png', 'names');
images_str = [images_str{:}];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create and prepare masks:
%__________________________________________________________________________
% make_mondrian_masks(sz_x, sz_y, n_masks, shape, selection)
% shape=2, circles. % selection=2, grayscale.
expmnt.num_masks =  ceil(length(frames_conds)/mask_framerate)*expmnt.n_trials_per_block.*2; % Minimum
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
if expmnt.debug_pt == 1
    WaitSecs(2);
else
    Screen('HideCursorHelper', w);
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

% if pf_estimation == 0

KbQueueCreate([], responseKeys);
KbQueueStart;
KbQueueFlush();
    trial_count = 0;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Start UML: before block begins
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pf_estimation == 1
        if expmnt.constant_stimuli_method == 1
            pf_stim_signals = expmnt.pf_stim_signals;
            trial_condition__pf_stim_signals = repmat(pf_stim_signals, 1, ceil(num_trials_block./length(pf_stim_signals)));
            if mod(num_reps_block, length(pf_stim_signals)) ~= 0
                % If number of repetitions of each trial type in each block
                % is not a multiple of the amount of pf_stim_signals defined
                % in settings, then, shuffle the orders, to have better chance
                % of having enough trials for each trial type:
                trial_condition__pf_stim_signals = Shuffle(trial_condition__pf_stim_signals);
            end
        else
            % Only create uml objects for conditions defined in 'expmnt.pf_type_estimations'
            simPhi0 = [0.60, 10, 0.5, 0.10]; % This line for debugging.

            if any(expmnt.pf_type_estimations == 1) % uml_pf.fam
                uml_pf.fam = UML(exp06_uml_settings('fam'));
                uml_pf.fam.setPhi0(simPhi0); % This line for debugging.
                uml_pf.fam.userdata01 = 'fam';
                pf_trial_count_fam = 0;
                pf_total_count_fam = 0;
            end

            if any(expmnt.pf_type_estimations == 2) % uml_pf.face
                uml_pf.face = UML(exp06_uml_settings('face'));
                uml_pf.face.setPhi0(simPhi0); % This line for debugging.
                uml_pf.face.userdata01 = 'face';
                pf_trial_count_face = 0;
                pf_total_count_face = 0;
            end

            if any(expmnt.pf_type_estimations == 3) % uml_pf.fam_abrupt
                uml_pf.fam_abrupt = UML(exp06_uml_settings('fam_abrupt'));
                uml_pf.fam_abrupt.setPhi0(simPhi0); % This line for debugging.
                uml_pf.fam_abrupt.userdata01 = 'fam_abrupt';
                pf_trial_count_fam_abrupt = 0;
                pf_total_count_fam_abrupt = 0;
            end
            
            if any(expmnt.pf_type_estimations == 4) % uml_pf.fam_faded
                uml_pf.fam_faded = UML(exp06_uml_settings('fam_faded'));
                uml_pf.fam_faded.setPhi0(simPhi0); % This line for debugging.
                uml_pf.fam_faded.userdata01 = 'fam_faded';
                pf_trial_count_fam_faded = 0;
                pf_total_count_fam_faded = 0;
            end

            % Codes in trial_condition__pf_type_estimations: 
            % 1 = fam; 2 = face PF curve estimation of the trial.
            % OR 3 = fam abrupt; 4 = fam faded onset 
            % Set the trial_condition__pf_type_estimations in the same order as the type of stimuli in "trial_condition__stimuli(1,:)".
            % ceil() and floor() are used just in case num of unfamfaces is odd.
%             trial_condition__pf_type_estimations = zeros(1, size(trial_condition__type, 2));
%             for type_of_trial = type_of_trial_index
%                 trial_condition__pf_type_estimations(trial_condition__type(1,:)==type_of_trial) = expmnt.pf_type_estimations(type_of_trial);
%             end
            clearvars type_of_trial;
        end
    end
    
for block = 1:length(expmnt.block_conds_order)

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Screen: stereogram to confirm vergence.
    %__________________________________________________________________________
    % Draw stereograms:
    Screen('DrawTextures', w, stereograms(which_stereogram(block+1)).tex_left, [], stereogram_box.left, []);
    Screen('DrawTextures', w, stereograms(which_stereogram(block+1)).tex_right, [], stereogram_box.right, []);
    Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
    Screen('TextSize', w, block_font_size); %
    DrawFormattedText(w, sprintf('Block: %i', block), screenXpixels*.01, screenYpixels-block_font_size, [1,1,1]); %posY_trial_info posX_trial_info
    Screen('Flip', w);
    WaitSecs(.5);
    if expmnt.savescreen; printscreenArray{end+1} = Screen('GetImage', w); end; %% Printscreen and save file.
    if expmnt.debug_pt == 1
        WaitSecs(3);
    else
        KbWait;
        WaitSecs(.1);
    end
    
    
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
    
    expmnt.trial_conds_order(block,:) =  Shuffle(expmnt.trial_conds_order(block,:));
    expmnt.identity_orders(block, :) =   Shuffle(expmnt.identity_orders(block,:));
    
    if pf_estimation
        % PF settings
        trial_condition__pf_type_estimations = Expand(expmnt.pf_type_estimations, expmnt.n_trials_per_block, 1);
        trial_condition__pf_type_estimations = Shuffle(trial_condition__pf_type_estimations);
        trial_condition__type =     [repmat(1, 1, length(unfamiliar_faces)),...
                                    repmat(2, 1, length(familiar_faces))];
        trial_index =               Shuffle(1:length(trial_condition__type));
    end
    
    %% Start trial frames:
    
    for trial = 1:expmnt.n_trials_per_block
        trial_count =           trial_count + 1;
        % Faces vector:
        frames_baseline = nan(1, ceil(length(frames_conds)/baseline_framerate));
        temp = [];
        idx_unfamiliar_faces = find(strcmp({images_str.cond}, 'unfamiliar') & strcmp({images_str.stim}, 'face'));
        idx_unfamiliar_houses = find(strcmp({images_str.cond}, 'unfamiliar') & strcmp({images_str.stim}, 'house'));
        idx_familiar_faces = find(strcmp({images_str.cond}, 'familiar') & strcmp({images_str.stim}, 'face'));
        
        % Every Xth image, present the oddball stimulus:
        idx_oddball = expmnt.odd_frequency:expmnt.odd_frequency:length(frames_baseline);
        
        switch expmnt.block_conds_order(block)
            case 0
                % PF estimation
                for i = 1:ceil(length(frames_baseline)/length([familiar_faces, unfamiliar_faces]))
                    temp = [temp, Shuffle([idx_unfamiliar_faces, idx_familiar_faces])];
                end
                % This is a vector with indices for each face to be presented:
                frames_baseline = temp(1:length(frames_baseline));
                clearvars temp;
                frames_baseline = frames_baseline(1);
                this_trial__images_cond = images_str(frames_baseline).cond;
                this_trial__stim_orientation = 0;
            case 1
                % Baseline: houses; oddball: unfamiliar faces
                for i = 1:ceil(length(frames_baseline)/length(unfamiliar_houses))
                    temp = [temp, Shuffle(idx_unfamiliar_houses)];
                end
                % This is a vector with indices for each face to be presented:
                frames_baseline = temp(1:length(frames_baseline));
                clearvars temp;
                temp = Shuffle(idx_unfamiliar_faces);
                frames_baseline(idx_oddball) = temp(1:length(frames_baseline(idx_oddball)));
                this_trial__stim_orientation = 0;
                this_trial__identity = 'none';
            case 2
                % Baseline: unfamiliar faces; oddball: familiar faces
                for i = 1:ceil(length(frames_baseline)/length(unfamiliar_faces))
                    temp = [temp, Shuffle(idx_unfamiliar_faces)];
                end
                % This is a vector with indices for each face to be presented:
                frames_baseline = temp(1:length(frames_baseline));
                clearvars temp;
                this_trial__identity = expmnt.identity_orders{block, trial};
                idx_correct_faces = find(strcmp({images_str.name}, this_trial__identity) & strcmp({images_str.cond}, 'familiar'));
                temp = [];
                for tmp_i = 1:(ceil(length(idx_oddball)/length(idx_correct_faces)))
                    temp = [temp, Shuffle(idx_correct_faces)];
                end
                frames_baseline(idx_oddball) = temp(1:length(frames_baseline(idx_oddball)));
                this_trial__stim_orientation = 0;
            case 3
                % Baseline: unfamiliar faces; oddball: familiar faces
                for i = 1:ceil(length(frames_baseline)/length(unfamiliar_faces))
                    temp = [temp, Shuffle(idx_unfamiliar_faces)];
                end
                % This is a vector with indices for each face to be presented:
                frames_baseline = temp(1:length(frames_baseline));
                clearvars temp;
                this_trial__identity = expmnt.identity_orders{block, trial};
                idx_correct_faces = find(strcmp({images_str.name}, this_trial__identity) & strcmp({images_str.cond}, 'familiar'));
                temp = [];
                for tmp_i = 1:(ceil(length(idx_oddball)/length(idx_correct_faces)))
                    temp = [temp, Shuffle(idx_correct_faces)];
                end
                frames_baseline(idx_oddball) = temp(1:length(frames_baseline(idx_oddball)));
                this_trial__stim_orientation = 180;
        end
        clearvars temp;

        % Expand each image index to the number of frames to be presented:
        frames_baseline = Expand(frames_baseline, baseline_framerate, 1);
        this_trial__images_idx = unique(frames_baseline, 'stable');
        this_trial__images_names = images_cell(this_trial__images_idx);

        % Define mask alpha: 0=visible, 1=suppressed
        if pf_estimation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % UML. Set masks alpha for this trial based on UML: begin trial
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If we are estimating PF, then overwrite alpha blending of masks:
            this_trial__pf_type_estimation       = trial_condition__pf_type_estimations(trial);
            if expmnt.constant_stimuli_method == 1
                this_trial__mask_alpha  = 1 - trial_condition__pf_stim_signals(trial);
            else
                if this_trial__pf_type_estimation == 1 
                    this_trial__uml_pf      = 'fam';
                    this_trial__mask_alpha  = 1 - uml_pf.fam.xnext;
                elseif  this_trial__pf_type_estimation == 2
                    this_trial__uml_pf      = 'face';
                    this_trial__mask_alpha  = 1 - uml_pf.face.xnext;
                elseif this_trial__pf_type_estimation == 3
                    this_trial__uml_pf      = 'fam_abrupt';
                    this_trial__mask_alpha  = 1 - uml_pf.fam_abrupt.xnext;
                elseif  this_trial__pf_type_estimation == 4
                    this_trial__uml_pf      = 'fam_faded';
                    this_trial__mask_alpha  = 1 - uml_pf.fam_faded.xnext;
                end
            end
        else
            switch expmnt.trial_conds_order(block,trial)
                case 0
                    this_trial__mask_alpha = expmnt.supp1__mask_alpha;
                case 1
                    this_trial__mask_alpha = expmnt.supp2__mask_alpha;
            end
        end
        
        
        % Mask position and potency:
        this_trial__xy_target_masks_displacement = ...
                    [expmnt.x_displacement, expmnt.y_displacement, expmnt.x_displacement, expmnt.y_displacement];
        this_trial__masks_position_xy =...
            mask_box_right + expmnt.cfs_1mirror_set_right_prop_of_displace * this_trial__xy_target_masks_displacement;

        % Target position:
        this_trial__target_position_xy = ...
            target_box_left - expmnt.cfs_1mirror_set_right_prop_of_displace .* this_trial__xy_target_masks_displacement.*[1,-1,1,-1];

        %% Trial loop
        %__________________________________________________________________________
        KbQueueFlush();

        % Construct Fade in:
        fade_in_frames = expmnt.fade_in_time * expmnt.monitor_hz;
        fade_in_cycle = round(linspace(0, this_trial__mask_alpha, fade_in_frames), 2);
        temp = Shuffle(1:expmnt.num_masks);
        temp = Expand(temp, mask_framerate, 1);
        fade_in_masks = temp(1:length(fade_in_cycle));
        clearvars temp;
        
        % Construct Fade out:
        fade_out_frames = expmnt.fade_out_time * expmnt.monitor_hz;
        fade_out_cycle = round(linspace(this_trial__mask_alpha, 0, fade_out_frames), 2);
        temp = Shuffle(1:expmnt.num_masks);
        temp = Expand(temp, mask_framerate, 1);
        fade_out_masks = temp(1:length(fade_out_cycle));
        
        % Run Fade in:
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
            Screen('DrawTextures', w, loadedStim(frames_baseline(frame)).theTexture, [],...
                this_trial__target_position_xy, this_trial__stim_orientation, [], 1-flipped(frame));
            % Draw masks:
            Screen('DrawTextures', w, masksTex(frames_masks(frame)).all, [],...
                this_trial__masks_position_xy, 180, [], this_trial__mask_alpha);
            vbl = Screen('Flip', w);
            timestamps(frame) = vbl;
        end
        % Run Fade out:
        for frame = 1:length(fade_out_cycle)
            Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
            % Draw masks:
            Screen('DrawTextures', w, masksTex(fade_out_masks(frame)).all, [],...
                this_trial__masks_position_xy, 180, [], fade_out_cycle(frame));
            Screen('Flip', w);
        end
        
        % Positions of identities in response screen:
        Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
        Screen('DrawTextures', w, loadedFam(1).theTexture, [],...
                responseItem_pos_all', this_trial__stim_orientation, [], 1);
        
%         Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
%         Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
%         Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
%         DrawFormattedText(w, 'Diaz Canel',  .50*screenXpixels, .25*screenYpixels, [1 1 1], [], [], [], [], [], textBox.right);
%         DrawFormattedText(w, 'Raul Castro', .75*screenXpixels, .50*screenYpixels, [1 1 1], [], [], [], [], [], textBox.right);
%         DrawFormattedText(w, 'Obama',       .50*screenXpixels, .75*screenYpixels, [1 1 1], [], [], [], [], [], textBox.right);
%         DrawFormattedText(w, 'Beyoncé',     .25*screenXpixels, .50*screenYpixels, [1 1 1], [], [], [], [], [], textBox.right);
        Screen('Flip', w);
        
        KbWait;
        
        % Gather response:
        [pressed, firstPress] = KbQueueCheck;
        if pressed
            if firstPress(keyUp)
                trial_response = 'DiazCanel';
                trial_response_time = firstPress(keyUp);
            end
            if firstPress(keyRight)
                trial_response = 'RaulCastro';
                trial_response_time = firstPress(keyRight);
            end
            if firstPress(keyDown)
                trial_response = 'BarackObama';
                trial_response_time = firstPress(keyDown);
            end
            if firstPress(keyLeft)
                trial_response = 'Beyonce';
                trial_response_time = firstPress(keyLeft);
            end
            if firstPress(spaceKey)
                response_correct = 0;
                trial_response = [];
                trial_response_time = [];
                break;
            end
            if firstPress(escapeKey)
                response_correct = 0;
                trial_response = [];
                trial_response_time = [];
                if firstPress(escapeKey)
                    closePT
                    return;
                end
                break;
            end
            response_correct = strcmp(trial_response, this_trial__identity);
            trial_responses = [trial_responses, trial_response];
            trial_response_times = [trial_response_times, trial_response_time];
        else
            response_correct = 0;
            trial_response_time = [];
        end
        
        if pf_estimation == 0
            WaitSecs(.2);
            Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
            Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
            Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
            DrawFormattedText(w, 'Confianza?',    0, .25, [1 1 1], [], [], [], [], [], textBox.right);
            DrawFormattedText(w, '1    2    3', 'center',   'center', [1 1 1], [], [], [], [], [], textBox.right);
            Screen('Flip', w);

            KbWait;

            % Gather response:
            [pressed, firstPress] = KbQueueCheck;
            if pressed
                if firstPress(keyLeft)
                    response_conf = 1;
                    response_conf_time = firstPress(keyLeft);
                end
                if firstPress(keyDown)
                    response_conf = 2;
                    response_conf_time = firstPress(keyDown);
                end
                if firstPress(keyRight)
                    response_conf = 3;
                    response_conf_time = firstPress(keyRight);
                end
            else
                response_conf = 0;
                response_conf_time = [];
            end
        end
        
        fprintf('\nTrial: %i, Mask alpha: %.03f, Correct: %i\n',...
            trial, this_trial__mask_alpha, response_correct);
        
    WaitSecs(.5);
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UML collection of response, updating function:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pf_estimation == 1
        % Keep track of valid trials:
        if expmnt.constant_stimuli_method == 1
            % Do nothing. Calculations are terminating the exp.
        else
            % Update UML object depending on PF function of this trial:
            if strcmp(this_trial__uml_pf, 'fam')
                    % response_correct_q2 = uml.simulateResponse(uml.xnext);
                    uml_pf.fam.update(response_correct);
                    fprintf('\n///////   PF estimation of this trial:    "fam"   ///////\n\n');
                    fprintf('\nResponse: %i | a = %.3f | b = %.3f | l = %.3f | xnext = %.3f',...
                        uml_pf.fam.r(end), uml_pf.fam.phi(end,1), uml_pf.fam.phi(end,2), uml_pf.fam.phi(end,4), uml_pf.fam.xnext);
                    fprintf('\nSweet points: %.3f | %.3f | %.3f | %.3f \n',...
                        uml_pf.fam.swpts(end, 1), uml_pf.fam.swpts(end, 2), uml_pf.fam.swpts(end, 3), uml_pf.fam.swpts(end, 4));
                if expmnt.debug_pt == 1
                    x = 0:.01:1;
                    figure; plot(x, ...
                    [0.5 + (1 - 0.5 - uml_pf.fam.phi(end,4)) .* (1 + exp(-(x - uml_pf.fam.phi(end,1)) .* uml_pf.fam.phi(end,2))).^(-1)]); % Current phi
                    xlim([0, 1]); ylim([0.50, 1]); title('Current estimation (phi)');
                end
                pf_trial_count_fam = pf_trial_count_fam + 1;
                pf_total_count_fam = pf_total_count_fam + 1;

            elseif strcmp(this_trial__uml_pf, 'face')
                    % response_correct_q1 = uml.simulateResponse(uml.xnext);
                    uml_pf.face.update(response_correct);
                    fprintf('\n///////   PF estimation of this trial:    "face"   ///////\n\n');
                    fprintf('\nResponse: %i | a = %.3f | b = %.3f | l = %.3f | xnext = %.3f',...
                        uml_pf.face.r(end), uml_pf.face.phi(end,1), uml_pf.face.phi(end,2), uml_pf.face.phi(end,4), uml_pf.face.xnext);
                    fprintf('\nSweet points: %.3f | %.3f | %.3f | %.3f \n',...
                        uml_pf.face.swpts(end, 1), uml_pf.face.swpts(end, 2), uml_pf.face.swpts(end, 3), uml_pf.face.swpts(end, 4));
                if expmnt.debug_pt == 1
                    x = 0:.01:1;
                    figure; plot(x, ...
                    [0.5 + (1 - 0.5 - uml_pf.face.phi(end,4)) .* (1 + exp(-(x - uml_pf.face.phi(end,1)) .* uml_pf.face.phi(end,2))).^(-1)]); % Current phi
                    xlim([0, 1]); ylim([0.50, 1]); title('Current estimation (phi)');
                end
                pf_trial_count_face = pf_trial_count_face + 1;
                pf_total_count_face = pf_total_count_face + 1;

            elseif strcmp(this_trial__uml_pf, 'fam_abrupt')
                    % response_correct_q2 = uml.simulateResponse(uml.xnext);
                    uml_pf.fam_abrupt.update(response_correct);
                    fprintf('\n///////   PF estimation of this trial:    "fam_abrupt"   ///////\n\n');
                    fprintf('\nResponse: %i | a = %.3f | b = %.3f | l = %.3f | xnext = %.3f',...
                        uml_pf.fam_abrupt.r(end), uml_pf.fam_abrupt.phi(end,1), uml_pf.fam_abrupt.phi(end,2), uml_pf.fam_abrupt.phi(end,4), uml_pf.fam_abrupt.xnext);
                    fprintf('\nSweet points: %.3f | %.3f | %.3f | %.3f \n',...
                        uml_pf.fam_abrupt.swpts(end, 1), uml_pf.fam_abrupt.swpts(end, 2), uml_pf.fam_abrupt.swpts(end, 3), uml_pf.fam_abrupt.swpts(end, 4));
                if expmnt.debug_pt == 1
                    x = 0:.01:1;
                    figure; plot(x, ...
                    [0.5 + (1 - 0.5 - uml_pf.fam_abrupt.phi(end,4)) .* (1 + exp(-(x - uml_pf.fam_abrupt.phi(end,1)) .* uml_pf.fam_abrupt.phi(end,2))).^(-1)]); % Current phi
                    xlim([0, 1]); ylim([0.50, 1]); title('Current estimation (phi)');
                end
                pf_trial_count_fam_abrupt = pf_trial_count_fam_abrupt + 1;
                pf_total_count_fam_abrupt = pf_total_count_fam_abrupt + 1;

            elseif strcmp(this_trial__uml_pf, 'fam_faded')
                    % response_correct_q2 = uml.simulateResponse(uml.xnext);
                    uml_pf.fam_faded.update(response_correct);
                    fprintf('\n///////   PF estimation of this trial:    "fam_faded"   ///////\n\n');
                    fprintf('\nResponse: %i | a = %.3f | b = %.3f | l = %.3f | xnext = %.3f',...
                        uml_pf.fam_faded.r(end), uml_pf.fam_faded.phi(end,1), uml_pf.fam_faded.phi(end,2), uml_pf.fam_faded.phi(end,4), uml_pf.fam_faded.xnext);
                    fprintf('\nSweet points: %.3f | %.3f | %.3f | %.3f \n',...
                        uml_pf.fam_faded.swpts(end, 1), uml_pf.fam_faded.swpts(end, 2), uml_pf.fam_faded.swpts(end, 3), uml_pf.fam_faded.swpts(end, 4));
                if expmnt.debug_pt == 1
                    x = 0:.01:1;
                    figure; plot(x, ...
                    [0.5 + (1 - 0.5 - uml_pf.fam_faded.phi(end,4)) .* (1 + exp(-(x - uml_pf.fam_faded.phi(end,1)) .* uml_pf.fam_faded.phi(end,2))).^(-1)]); % Current phi
                    xlim([0, 1]); ylim([0.50, 1]); title('Current estimation (phi)');
                end
                pf_trial_count_fam_faded = pf_trial_count_fam_faded + 1;
                pf_total_count_fam_faded = pf_total_count_fam_faded + 1;

            end
        end
    end

    %% Save data
    %__________________________________________________________________________
    data(trial_count).subject =               subject;
    data(trial_count).block =                 block;
    data(trial_count).trial =                 trial;
    data(trial_count).images_idx =            this_trial__images_idx;
    data(trial_count).images_idx =            this_trial__images_names;
    data(trial_count).timestamps =            timestamps;
    data(trial_count).stim_orientation =      this_trial__stim_orientation;
    data(trial_count).this_trial__identity =  this_trial__identity;
    data(trial_count).trial_response =        trial_response;
    data(trial_count).trial_response_time =   trial_response_time;
    data(trial_count).response_correct =      response_correct;
    if pf_estimation == 0
        data(trial_count).response_conf =         response_conf;
        data(trial_count).response_conf_time =    response_conf_time;
    end
    data(trial_count).this_trial__mask_alpha = this_trial__mask_alpha;
    data(trial_count).trial_responses =       trial_responses;
    data(trial_count).trial_response_times =  trial_response_times;
    data(trial_count).frames_conds =          frames_conds;
    data(trial_count).frames_baseline =       frames_baseline;
    data(trial_count).frames_masks =          frames_masks;
    
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

    
    %% Save pf estimation
    if ~exist('data', 'var')
        data = [];
    end
    if pf_estimation == 1
        save([filename, '_beh.mat'], 'data', 'expmnt', 'uml_pf');
    else
        save([filename, '_beh.mat'], 'data', 'expmnt');
    end
    
end
    
% end
    
%% _________________________________________________________________________
%|
%|
%|
%|
%|                    BEGINNING OF PF 
%|
%|
%|
%|__________________________________________________________________________

if 0 % pf_estimation
% PF settings
expmnt.constant_stimuli_method = 0;
expmnt.pf_type_estimations = 1;
expmnt.pf_trials_each_per_block = 35;
expmnt.block_pf = 3;

for block = 1:expmnt.block_pf

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
    frames_baseline = nan(1, length(frames_conds)/baseline_framerate);
    temp = [];
    for i = 1:ceil(length(frames_baseline)/length(unfamiliar_faces))
        temp = [temp, Shuffle(1:length(unfamiliar_faces))];
    end
    % This is a vector with indices for each face to be presented:
    frames_baseline = temp(1:length(frames_baseline));
    clearvars temp;
    % Every Xth image, present the oddball stimulus:
    idx_oddball = expmnt.odd_frequency:expmnt.odd_frequency:length(frames_baseline);
    temp = Shuffle((1:length(familiar_faces))+length(unfamiliar_faces));
    frames_baseline(idx_oddball) = temp(1:length(frames_baseline(idx_oddball)));
    clearvars temp;

    % Expand each image index to the number of frames to be presented:
    frames_baseline = Expand(frames_baseline, baseline_framerate, 1);
    
    % Define mask alpha: 0=visible, 1=suppressed
    switch expmnt.block_conds_order(block)
        case 0
            this_trial__mask_alpha = expmnt.supp1__mask_alpha;
        case 1
            this_trial__mask_alpha = expmnt.supp2__mask_alpha;
    end
    
    % Mask position and potency:
    this_trial__xy_target_masks_displacement = ...
                [expmnt.x_displacement, expmnt.y_displacement, expmnt.x_displacement, expmnt.y_displacement];
    this_trial__masks_position_xy =...
        mask_box_right + expmnt.cfs_1mirror_set_right_prop_of_displace * this_trial__xy_target_masks_displacement;

    % Target position:
    this_trial__target_position_xy = ...
        target_box_left - expmnt.cfs_1mirror_set_right_prop_of_displace .* this_trial__xy_target_masks_displacement.*[1,-1,1,-1];
    
    % PF settings
    trial_condition__pf_type_estimations = Expand(expmnt.pf_type_estimations, expmnt.pf_trials_each_per_block, 1);
    trial_condition__pf_type_estimations = Shuffle(trial_condition__pf_type_estimations);
    trial_condition__type =     [repmat(1, 1, length(unfamiliar_faces)),...
                                repmat(2, 1, length(familiar_faces))];
    trial_index =               Shuffle(1:length(trial_condition__type));

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Start UML: before block begins
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pf_estimation == 1
        if expmnt.constant_stimuli_method == 1
            pf_stim_signals = expmnt.pf_stim_signals;
            trial_condition__pf_stim_signals = repmat(pf_stim_signals, 1, ceil(num_trials_block./length(pf_stim_signals)));
            if mod(num_reps_block, length(pf_stim_signals)) ~= 0
                % If number of repetitions of each trial type in each block
                % is not a multiple of the amount of pf_stim_signals defined
                % in settings, then, shuffle the orders, to have better chance
                % of having enough trials for each trial type:
                trial_condition__pf_stim_signals = Shuffle(trial_condition__pf_stim_signals);
            end
        else
            % Only create uml objects for conditions defined in 'expmnt.pf_type_estimations'
            simPhi0 = [0.60, 10, 0.5, 0.10]; % This line for debugging.

            if any(expmnt.pf_type_estimations == 1) % uml_pf.fam
                uml_pf.fam = UML(exp06_uml_settings('fam'));
                uml_pf.fam.setPhi0(simPhi0); % This line for debugging.
                uml_pf.fam.userdata01 = 'fam';
                pf_trial_count_fam = 0;
                pf_total_count_fam = 0;
            end

            if any(expmnt.pf_type_estimations == 2) % uml_pf.face
                uml_pf.face = UML(exp06_uml_settings('face'));
                uml_pf.face.setPhi0(simPhi0); % This line for debugging.
                uml_pf.face.userdata01 = 'face';
                pf_trial_count_face = 0;
                pf_total_count_face = 0;
            end

            if any(expmnt.pf_type_estimations == 3) % uml_pf.fam_abrupt
                uml_pf.fam_abrupt = UML(exp06_uml_settings('fam_abrupt'));
                uml_pf.fam_abrupt.setPhi0(simPhi0); % This line for debugging.
                uml_pf.fam_abrupt.userdata01 = 'fam_abrupt';
                pf_trial_count_fam_abrupt = 0;
                pf_total_count_fam_abrupt = 0;
            end
            
            if any(expmnt.pf_type_estimations == 4) % uml_pf.fam_faded
                uml_pf.fam_faded = UML(exp06_uml_settings('fam_faded'));
                uml_pf.fam_faded.setPhi0(simPhi0); % This line for debugging.
                uml_pf.fam_faded.userdata01 = 'fam_faded';
                pf_trial_count_fam_faded = 0;
                pf_total_count_fam_faded = 0;
            end

            % Codes in trial_condition__pf_type_estimations: 
            % 1 = fam; 2 = face PF curve estimation of the trial.
            % OR 3 = fam abrupt; 4 = fam faded onset 
            % Set the trial_condition__pf_type_estimations in the same order as the type of stimuli in "trial_condition__stimuli(1,:)".
            % ceil() and floor() are used just in case num of unfamfaces is odd.
%             trial_condition__pf_type_estimations = zeros(1, size(trial_condition__type, 2));
%             for type_of_trial = type_of_trial_index
%                 trial_condition__pf_type_estimations(trial_condition__type(1,:)==type_of_trial) = expmnt.pf_type_estimations(type_of_trial);
%             end
            clearvars type_of_trial;
        end
    end
    
    
    %% select alpha for the trial
    for this_trial__id = 1:length(trial_condition__pf_type_estimations)

        if pf_estimation == 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % UML. Set masks alpha for this trial based on UML: begin trial
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If we are estimating PF, then overwrite alpha blending of masks:
            this_trial__pf_type_estimation       = trial_condition__pf_type_estimations(this_trial__id);
            if expmnt.constant_stimuli_method == 1
                this_trial__mask_alpha  = 1 - trial_condition__pf_stim_signals(this_trial__id);
            else
                if this_trial__pf_type_estimation == 1 
                    this_trial__uml_pf      = 'fam';
                    this_trial__mask_alpha  = 1 - uml_pf.fam.xnext;
                elseif  this_trial__pf_type_estimation == 2
                    this_trial__uml_pf      = 'face';
                    this_trial__mask_alpha  = 1 - uml_pf.face.xnext;
                elseif this_trial__pf_type_estimation == 3
                    this_trial__uml_pf      = 'fam_abrupt';
                    this_trial__mask_alpha  = 1 - uml_pf.fam_abrupt.xnext;
                elseif  this_trial__pf_type_estimation == 4
                    this_trial__uml_pf      = 'fam_faded';
                    this_trial__mask_alpha  = 1 - uml_pf.fam_faded.xnext;
                end
            end
        end


        %% Trial
        expmnt.fade_in_time = 0.5;
        expmnt.fade_out_time = 0.5;
        stim_for_this_trial = trial_index(this_trial__id);
        stim_condition_for_this_trial = trial_condition__type(trial_index(this_trial__id));

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

        for frame = 1:baseline_framerate
            Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
            % Draw target:
            Screen('DrawTextures', w, loadedStim(stim_for_this_trial).theTexture, [],...
                this_trial__target_position_xy, 0, [], 1-flipped(frame));
            % Draw masks:
            Screen('DrawTextures', w, masksTex(frames_masks(frame)).all, [],...
                this_trial__masks_position_xy, 180, [], this_trial__mask_alpha);

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
        
        Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
        Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
        Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
        Screen('Flip', w);
        WaitSecs(1);

        % Gather response:
        [pressed, firstPress] = KbQueueCheck;
        if pressed
            if firstPress(key1)
                trial_response = 1;
                trial_response_time = firstPress(key1);
            end
            if firstPress(key2)
                trial_response = 2;
                trial_response_time = firstPress(key2);
            end
            if firstPress(spaceKey)
                break;
            end
            trial_responses = [trial_responses, trial_response];
            trial_response_times = [trial_response_times, trial_response_time];
            response_correct_q2 = stim_condition_for_this_trial == trial_response;
        else
            response_correct_q2 = 0;
        end

        fprintf('\nTrial: %i, Estimating: %i, Mask alpha: %.03f\n',...
            this_trial__id, this_trial__pf_type_estimation, this_trial__mask_alpha);
        
        if response_correct_q2
            Screen('DrawTextures', w, textsImages(6).textures, [], [textBox.left; textBox.right]', 0, [], []);
            Screen('DrawLines', w, fixCoordsLeft, fix_lineWidth, fix_color, fix_position_left, 2);
            Screen('DrawLines', w, fixCoordsRight, fix_lineWidth, fix_color, fix_position_right, 2);
            Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
            Screen('Flip', w);
            WaitSecs(.5);
        end


        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UML collection of response, updating function:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pf_estimation == 1
            % Keep track of valid trials:
            if expmnt.constant_stimuli_method == 1
                % Do nothing. Calculations are terminating the exp.
            else
                % Update UML object depending on PF function of this trial:
                if strcmp(this_trial__uml_pf, 'fam')
                        % response_correct_q2 = uml.simulateResponse(uml.xnext);
                        uml_pf.fam.update(response_correct_q2);
                        fprintf('\n///////   PF estimation of this trial:    "fam"   ///////\n\n');
                        fprintf('\nResponse: %i | a = %.3f | b = %.3f | l = %.3f | xnext = %.3f',...
                            uml_pf.fam.r(end), uml_pf.fam.phi(end,1), uml_pf.fam.phi(end,2), uml_pf.fam.phi(end,4), uml_pf.fam.xnext);
                        fprintf('\nSweet points: %.3f | %.3f | %.3f | %.3f \n',...
                            uml_pf.fam.swpts(end, 1), uml_pf.fam.swpts(end, 2), uml_pf.fam.swpts(end, 3), uml_pf.fam.swpts(end, 4));
                    if expmnt.debug_pt == 1
                        x = 0:.01:1;
                        figure; plot(x, ...
                        [0.5 + (1 - 0.5 - uml_pf.fam.phi(end,4)) .* (1 + exp(-(x - uml_pf.fam.phi(end,1)) .* uml_pf.fam.phi(end,2))).^(-1)]); % Current phi
                        xlim([0, 1]); ylim([0.50, 1]); title('Current estimation (phi)');
                    end
                    pf_trial_count_fam = pf_trial_count_fam + 1;
                    pf_total_count_fam = pf_total_count_fam + 1;

                elseif strcmp(this_trial__uml_pf, 'face')
                        % response_correct_q1 = uml.simulateResponse(uml.xnext);
                        uml_pf.face.update(response_correct_q1);
                        fprintf('\n///////   PF estimation of this trial:    "face"   ///////\n\n');
                        fprintf('\nResponse: %i | a = %.3f | b = %.3f | l = %.3f | xnext = %.3f',...
                            uml_pf.face.r(end), uml_pf.face.phi(end,1), uml_pf.face.phi(end,2), uml_pf.face.phi(end,4), uml_pf.face.xnext);
                        fprintf('\nSweet points: %.3f | %.3f | %.3f | %.3f \n',...
                            uml_pf.face.swpts(end, 1), uml_pf.face.swpts(end, 2), uml_pf.face.swpts(end, 3), uml_pf.face.swpts(end, 4));
                    if expmnt.debug_pt == 1
                        x = 0:.01:1;
                        figure; plot(x, ...
                        [0.5 + (1 - 0.5 - uml_pf.face.phi(end,4)) .* (1 + exp(-(x - uml_pf.face.phi(end,1)) .* uml_pf.face.phi(end,2))).^(-1)]); % Current phi
                        xlim([0, 1]); ylim([0.50, 1]); title('Current estimation (phi)');
                    end
                    pf_trial_count_face = pf_trial_count_face + 1;
                    pf_total_count_face = pf_total_count_face + 1;

                elseif strcmp(this_trial__uml_pf, 'fam_abrupt')
                        % response_correct_q2 = uml.simulateResponse(uml.xnext);
                        uml_pf.fam_abrupt.update(response_correct_q2);
                        fprintf('\n///////   PF estimation of this trial:    "fam_abrupt"   ///////\n\n');
                        fprintf('\nResponse: %i | a = %.3f | b = %.3f | l = %.3f | xnext = %.3f',...
                            uml_pf.fam_abrupt.r(end), uml_pf.fam_abrupt.phi(end,1), uml_pf.fam_abrupt.phi(end,2), uml_pf.fam_abrupt.phi(end,4), uml_pf.fam_abrupt.xnext);
                        fprintf('\nSweet points: %.3f | %.3f | %.3f | %.3f \n',...
                            uml_pf.fam_abrupt.swpts(end, 1), uml_pf.fam_abrupt.swpts(end, 2), uml_pf.fam_abrupt.swpts(end, 3), uml_pf.fam_abrupt.swpts(end, 4));
                    if expmnt.debug_pt == 1
                        x = 0:.01:1;
                        figure; plot(x, ...
                        [0.5 + (1 - 0.5 - uml_pf.fam_abrupt.phi(end,4)) .* (1 + exp(-(x - uml_pf.fam_abrupt.phi(end,1)) .* uml_pf.fam_abrupt.phi(end,2))).^(-1)]); % Current phi
                        xlim([0, 1]); ylim([0.50, 1]); title('Current estimation (phi)');
                    end
                    pf_trial_count_fam_abrupt = pf_trial_count_fam_abrupt + 1;
                    pf_total_count_fam_abrupt = pf_total_count_fam_abrupt + 1;

                elseif strcmp(this_trial__uml_pf, 'fam_faded')
                        % response_correct_q2 = uml.simulateResponse(uml.xnext);
                        uml_pf.fam_faded.update(response_correct_q2);
                        fprintf('\n///////   PF estimation of this trial:    "fam_faded"   ///////\n\n');
                        fprintf('\nResponse: %i | a = %.3f | b = %.3f | l = %.3f | xnext = %.3f',...
                            uml_pf.fam_faded.r(end), uml_pf.fam_faded.phi(end,1), uml_pf.fam_faded.phi(end,2), uml_pf.fam_faded.phi(end,4), uml_pf.fam_faded.xnext);
                        fprintf('\nSweet points: %.3f | %.3f | %.3f | %.3f \n',...
                            uml_pf.fam_faded.swpts(end, 1), uml_pf.fam_faded.swpts(end, 2), uml_pf.fam_faded.swpts(end, 3), uml_pf.fam_faded.swpts(end, 4));
                    if expmnt.debug_pt == 1
                        x = 0:.01:1;
                        figure; plot(x, ...
                        [0.5 + (1 - 0.5 - uml_pf.fam_faded.phi(end,4)) .* (1 + exp(-(x - uml_pf.fam_faded.phi(end,1)) .* uml_pf.fam_faded.phi(end,2))).^(-1)]); % Current phi
                        xlim([0, 1]); ylim([0.50, 1]); title('Current estimation (phi)');
                    end
                    pf_trial_count_fam_faded = pf_trial_count_fam_faded + 1;
                    pf_total_count_fam_faded = pf_total_count_fam_faded + 1;

                end
            end
        end

        %% Save pf estimation
        if pf_estimation == 1
            save([filename, '_beh.mat'], 'uml_pf');
        end

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
    
end

    
end

%% Close PsychToolbox
%__________________________________________________________________________

closePT;

%% Plot test of cycles:
%__________________________________________________________________________
if pf_estimation == 0
    figure;
    subplot(4,1,1);
    plot(data(1).frames_baseline<73); 
    title(sprintf('Unfamiliar faces: %.02f Hz', expmnt.baseline_hz));
    ylim([-0.1, 1.1]);
    subplot(4,1,2);
    plot(data(1).frames_baseline>72); 
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
    % plot(data.frames_baseline);
    % title('Unfamiliar faces: 1-72, Familiar faces: 73-144');

end

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
    fix_position_left =  [xCenter - x_displacement_left,  yCenter + expmnt.y_displacement]; 
    fix_position_right = [xCenter + x_displacement_right, yCenter + expmnt.y_displacement]; 
    
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
        vergence_bar_left =     vergence_bar_left - [x_displacement_left, -expmnt.y_displacement,...
            x_displacement_left, -expmnt.y_displacement];
        vergence_bar_right =    vergence_bar_right + [x_displacement_right, expmnt.y_displacement,...
            x_displacement_right, expmnt.y_displacement];
        vergence_bars_positions = [ vergence_bar_left + [vergence_bar_offplace_left, 0, vergence_bar_offplace_left, 0];...
                                    vergence_bar_left - [vergence_bar_offplace_left, 0, vergence_bar_offplace_left, 0];...
                                    vergence_bar_right + [vergence_bar_offplace_right, 0, vergence_bar_offplace_right, 0];...
                                    vergence_bar_right - [vergence_bar_offplace_right, 0, vergence_bar_offplace_right, 0] ]';
    else
        vergence_bar_offplace = expmnt.mask_width/2 + (expmnt.vergence_bar_width*1.25/2);
		vergence_bar = CenterRectOnPointd(expmnt.vergence_bar, xCenter, yCenter);
		vergence_bar_right = vergence_bar + [expmnt.x_displacement, expmnt.y_displacement,...
            expmnt.x_displacement, expmnt.y_displacement];
		vergence_bar_left = vergence_bar - [expmnt.x_displacement, -expmnt.y_displacement,...
            expmnt.x_displacement, -expmnt.y_displacement];
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
	textBox.right = textBox_right + [x_displacement_right, expmnt.y_displacement, x_displacement_right, expmnt.y_displacement];
	textBox.left =  textBox_left - [x_displacement_left, -expmnt.y_displacement, x_displacement_left, -expmnt.y_displacement];
    
    
    % =========================================================================
    % Stereograms:
    % =========================================================================
    stereogram_box_left = CenterRectOnPointd(mask_box*expmnt.cfs_1mirror_set_size_ratio, xCenter, yCenter); 
    stereogram_box_right = CenterRectOnPointd(mask_box, xCenter, yCenter); 
    stereogram_box.left = stereogram_box_left - [x_displacement_left,  -expmnt.y_displacement,...
        x_displacement_left, -expmnt.y_displacement];
    stereogram_box.right = stereogram_box_right  + [x_displacement_right, expmnt.y_displacement,...
        x_displacement_right, expmnt.y_displacement];
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

%% Calculate alpha of masks from PF:
function [masks_alpha_blending] = calculate_alpha_blending_threshold(plot_pf_function, uml_object, perf_over_gamma)
    % calculate_alpha_blending_threshold(expmnt, matPFFile)
    %
    % p_log = ( ( (1-lambda) - gamma ) * y ) + gamma;
    p_log_for_supp  = ( ( (1-uml_object.phi(end,4)) - uml_object.phi(end,3) ) * perf_over_gamma ) + uml_object.phi(end,3);
    signal_x        = -log(((1 - uml_object.phi(end,3) - uml_object.phi(end,4)) / (p_log_for_supp - uml_object.phi(end,3))) - 1)*(1/uml_object.phi(end,2)) + uml_object.phi(end,1);
    % Change direction: from 'amount of signal' to 'mask alpha blending'
    masks_alpha_blending = 1 - signal_x;
    % Plot for debugging:
    if plot_pf_function == 1
        x_values = 0:0.01:1;
        plot_fam = plot(x_values, uml_object.phi(end,3)+(1-uml_object.phi(end,3)-uml_object.phi(end,4)).*(1+exp(-(x_values-uml_object.phi(end,1)).*uml_object.phi(end,2))).^(-1)); % Current phi
        line([signal_x, signal_x], [0.50, 1]);
        xlim([0, 1]); ylim([0.50, 1]);
    end
    fprintf('\n=====================================================');
    fprintf('\nPF %s: alpha= %.2f | beta= %.2f | lambda= %.2f |', uml_object.userdata01, uml_object.phi(end,1), uml_object.phi(end,2), uml_object.phi(end,4));
end

%% Function to close PsychToolbox cleanly:
function [] = closePT()

    % Release keyboard with PsychToolbox, wait and clear workspace:
    KbQueueRelease();
    KbReleaseWait();
    WaitSecs(1);
    sca;
    clear Screen;
    ListenChar(0); % Restore keyboard output to Matlab

end

%% Done:
%__________________________________________________________________________
% March 2024:
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
% [X] Adjust mask and faces sizes: based on visual angle degrees.
% [X] Add selection of short trials: 10 famfaces cycles? Change in duration
% [X] Code UML psychophysics estimation.
% [X] Set parts mask alpha from UML estimation.
%  X  Test file pf_estimation_block
%  X  Then insert PsychToolbox trial presentation.
%  X  Then insert the code into this script.
% June 2024:
% [X] Can we use Webmorph for automatic face images processing?
%  X  Yes, we can. At least with rigid align.
% [-] Add task: categorize Cneuro vs famous?
%     Define a specific face category (females), specify positions in the
%     block, and score if: cneuro vs famous.
% [X] Later: clean ways that PF and main task use for loop.
%     Main task uses the frames_cond, but PF uses trial level.
% [X] Set custom high.
% [X] What can work as *divider*?
%     Solution: Use the tall table, and construct a cardboard box with two
%     tunnels (each for each eye), that can be adjusted in tilt and height.
%     Important that the material does not reflect too much.

