
%% Get features of each face from filenames:

images_list = dir('C:\Users\cthul\Dropbox\PROYECTOS\0-PhD\ssVEP_CFS\stimuli\caras_experimento\*.png')
images_list = {images_list.name}

regexp('familiar_cneuro_m_Adrian.png',...
    '(?<cond>\w+)_(?<where>\w+)_(?<sex>\w+)_(?<name>\w+).png', 'names')

images_str = regexp(images_list,...
    '(?<cond>\w+)_(?<where>\w+)_(?<sex>\w+)_(?<name>\w+).png', 'names')
images_str = [images_str{:}];
images_str(1).cond
images_str(1).where
images_str(1).sex
images_str(1).name


%%
