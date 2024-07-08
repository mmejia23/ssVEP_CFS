
%% Debug: -----------------------------------------------------------------
% Usar estas dos lineas para 'probar' el codigo. Desactiva el trigger al
% BrainVission, y muestra la pantalla transparente.
% Threshold estimation
[data, expmnt] = ssvepCFS('NombreParticipante', 'pf',1, 'debug',1, 'triggers',0);
% Main task
[data, expmnt] = ssvepCFS('NombreParticipante', 'debug',1, 'triggers',0);

%% Run: -------------------------------------------------------------------
% 1-Primero se estima el umbral para el participante ('pf', 1).
% 2-Segundo, se ejecuta la tarea principal completa.
%% Task 1: Psychophysics task to estimate alpha blending values:
% Stereograms En orden: 1=A; 2=C; 3=H; 4=J;
[data, expmnt] = ssvepCFS('NombreParticipante1', 'pf',1);

%% Task 2: After estimation of alpha blending values:
% Stereograms En orden: 1=L; 2=O; 3=P; 4=S; 5=T; 6=X; 7=Z;
ssvepCFS('NombreParticipante', 'calculate_alpha');
[data, expmnt] = ssvepCFS('NombreParticipante1', 'blockstart',1);

