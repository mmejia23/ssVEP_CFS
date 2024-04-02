


% [data, expmnt] = ssvepCFS(subject, suppression, supp_alpha);
% suppression: 0=visible,  1=suprimido
% supp_alpha:  anotar numero de curva psicofisica

% Psychophysics task to estimate alpha blending values:


% Knowing alpha blending values:
[data, expmnt] = ssvepCFS('ManuelMejia', [0.04, 0.32])

% Not knowing alpha blending values:
[data, expmnt] = ssvepCFS('ManuelMejia')

