%% Calculate signal value that gives us a specific performance:
function [masks_alpha_blending] = calculate_signal_threshold(uml_object, perf)
    % Invert function
%     uml_pf.fam.psycfun

    x_min = 0.03;
    y_max = uml_object.psycfun(1-x_min,...
        uml_object.phi(end,1), uml_object.phi(end,2),...
        uml_object.phi(end,3), uml_object.phi(end,4));

    % Set target value on Y:
%     perf = 0.53
    % Set a guess for the value of X, that corresponds to Y:
    guess = 0.5;
    % Invert the function, with defined values on the other variables:a,b,g,l
    % The 'inversion' is done by the trick of having the function to be zero
    % for the value that we want: function result and gval are equal
    inv_fun = @(x) (uml_object.psycfun(x,...
        uml_object.phi(end,1), uml_object.phi(end,2),...
        uml_object.phi(end,3), uml_object.phi(end,4)) - perf);
    % Run fzero that tries to find the x value that gives you zero
    signal_x = fzero(@(x) inv_fun(x), guess);
    inv_fun(signal_x);
    % If you use the resulting S value, then you can confirm that it gives you
    % the Y value that you wanted:
    uml_object.psycfun(signal_x,...
        uml_object.phi(end,1), uml_object.phi(end,2),...
        uml_object.phi(end,3), uml_object.phi(end,4));
    masks_alpha_blending = 1 - signal_x;
end
