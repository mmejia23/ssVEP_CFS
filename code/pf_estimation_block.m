pf_estimation = 1;
expmnt.debug_pt = 0;
expmnt.constant_stimuli_method = 0;
expmnt.pf_type_estimations = 1;
expmnt.pf_trials_each = 80;
trial_condition__pf_type_estimations = Expand(expmnt.pf_type_estimations, expmnt.pf_trials_each, 1);
trial_condition__pf_type_estimations = Shuffle(trial_condition__pf_type_estimations);
trial_condition__type =     [repmat(1, 1, length(unfamiliar_images)),...
                            repmat(2, 1, length(familiar_images))];
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

    for frame = 1:length(frames_conds)

        Screen('DrawTextures', w, vergence(1).tex, [], vergence_bars_positions, 0, [], []);
        % Draw target:
        Screen('DrawTextures', w, loadedFaces(stim_for_this_trial).theTexture, [],...
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
        save([filename, 'session-pf_beh.mat'], 'uml_pf');
    end
    
end
    