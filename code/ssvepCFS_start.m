
%% Debug: -----------------------------------------------------------------
% Usar estas dos lineas para 'probar' el codigo. Desactiva el trigger al
% BrainVission, y muestra la pantalla transparente.
% Estimation of threshold
[data, expmnt] = ssvepCFS('NombreParticipante', 'pf', 1, 'debug',1, 'triggers',0);

% Main task
[data, expmnt] = ssvepCFS('NombreParticipante', 'alphas', [0.03,0.30], 'debug',1, 'triggers',0);



%% Run: -------------------------------------------------------------------
% Usar estos codigos para ejecutar la tarea.
% 1-Primero se estima el umbral para el participante ('pf', 1).
% 2-Luego se calculan los valores para el alpha blending 'calculate_alpha'
% 3-Y por ultimo, se ejecuta la tarea principal completa: pero agregar los
% valores que den para el alpha (o usar los default para probar).
%% Task 1: Psychophysics task to estimate alpha blending values:
% Stereograms En orden: 1=A; 2=C; 3=H; 4=P; 5=L; 6=S; 7=T; 8=X; 9=Z;
[data, expmnt] = ssvepCFS('NombreParticipante', 'pf', 1, 'sim', 1);
 
%% Task 2: Knowing alpha blending values:
% Stereograms En orden: 1=A; 2=C; 3=H; 4=P; 5=L; 6=S; 7=T; 8=X; 9=Z;
% 2=Upright, 3=Inverted
ssvepCFS('NombreParticipante', 'calculate_alpha');
[data, expmnt] = ssvepCFS('NombreParticipante', 'blockstart', 1, 'sim', 1);



