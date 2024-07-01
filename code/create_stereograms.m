addpath(genpath('C:\Users\cthul\Dropbox\PROYECTOS\0-PhD\MLtoolboxes\YkImageTools'));

temp_img = imread("C:\Users\cthul\Downloads\pngtree-by-black-and-white-sign-png-image_4583570.png");
temp_img = double(temp_img);
temp_img(temp_img<255/2) = 15;
temp_img(temp_img>255/2) = 0;
[im_left, im_right] = stereogram2(temp_img);
imwrite(im_left, "C:\Users\cthul\Downloads\stereogram-01_left.png");
imwrite(im_right, "C:\Users\cthul\Downloads\stereogram-01_right.png");

dir_to_save = 'C:\Users\cthul\Dropbox\PROYECTOS\0-PhD\ssVEP_CFS\stimuli\stereograms2';
stereogram_letters = {'A', 'C', 'H', 'J', 'L', 'M', 'O', 'P', 'S', 'T', 'U', 'W', 'X', 'Y', 'Z'};
for i = 1:length(stereogram_letters)
    letter = stereogram_letters{i};
    temp_img = insertText(zeros(150,150),[75 75],letter,...
        'BoxOpacity',0,...
        'FontSize',140,'TextColor','w', 'AnchorPoint', 'center',...
        'Font', 'LucidaTypewriterBold');
    temp_img = sum(temp_img, 3);
    temp_img(temp_img<max(temp_img, [], 'all')/2) = 15;
    temp_img(temp_img>max(temp_img, [], 'all')/2) = 0;
    figure; imshow(temp_img);
    [im_left, im_right] = stereogram2(temp_img, 'plot');
    imwrite(im_left, fullfile(dir_to_save, sprintf('stereogram_%02d_%s_l.png', i, letter)));
    imwrite(im_right, fullfile(dir_to_save, sprintf('stereogram_%02d_%s_r.png', i, letter)));
end

