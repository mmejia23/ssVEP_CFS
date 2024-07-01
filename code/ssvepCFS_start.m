


% [data, expmnt] = ssvepCFS(subject, suppression, supp_alpha);
% suppression: 0=visible,  1=suprimido
% supp_alpha:  anotar numero de curva psicofisica

%% Psychophysics task to estimate alpha blending values:
[data, expmnt] = ssvepCFS('ManuelMejia', 'pf', 1) 
 
% Knowing alpha blending values:
 
%% Not knowing alpha blending values:

[data, expmnt] = ssvepCFS('ManuelMejia', 'alphas', [0.03, 0.3], 'debug',0)
[data, expmnt] = ssvepCFS('ManuelMejia', 'alphas', [0.03, 0.3])



