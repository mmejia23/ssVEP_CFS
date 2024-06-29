


% [data, expmnt] = ssvepCFS(subject, suppression, supp_alpha);
% suppression: 0=visible,  1=suprimido
% supp_alpha:  anotar numero de curva psicofisica

%% Psychophysics task to estimate alpha blending values:
[data, expmnt] = ssvepCFS('ManuelMejia', 'pf', 1) 
 
% Knowing alpha blending values:
 
%% Not knowing alpha blending values:

[data, expmnt] = ssvepCFS('ManuelMejia', 'alphas', [0.03, 0.3], 'debug',1)
[data, expmnt] = ssvepCFS('ManuelMejia', 'alphas', [0.03, 0.3])


    
%% 
%{
[ ] Fotos completas
[ ] Condicion de casas
[ ] Que la ultima no sea familiar
[ ] Agregar la cantidad de estimulos para igualar a Rossion
[ ] Programar condicion con un solo estimulo central

[ ] Programar condicion sin CFS

%}