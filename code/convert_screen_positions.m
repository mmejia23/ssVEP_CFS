function [result, center] = convert_screen_positions(direction, x, image, n_images, screen)
% Function to convert between an index of a list of images presented on
% screen, to their box coordinates for Psychtoolbox. And viceversa: from
% the coordinates clicked on screen, to the index of the image clicked.
if direction == 0
    n_col = round(n_images ./ (screen(1)./ screen(2)));
    n_row = ceil(n_images ./ n_col);
    area_for_each_image = (screen(1).*screen(2)) ./ (n_row.*n_col);
    if image(1)<image(2)
        img_ratio = image(2)/image(1);
        adj_width = floor(sqrt(area_for_each_image ./ img_ratio));
        adj_height = floor(adj_width.*img_ratio);
    else
        img_ratio = image(2)/image(1);
        adj_width = floor(screen(1) ./ n_col);
        adj_height = floor(adj_width.*img_ratio);
    end
    img_index = x;
%     rows(1) = floor(screen(1)./adj_width);
    x_row = rem(img_index-1, n_col) + 1;
    y_row = ceil(img_index /n_col);
    % Get the box for each image:
    result = [ 0, 0, adj_width, adj_height ];
    center = round([ x_row*adj_width-adj_width./2, y_row*adj_height-adj_height./2 ]);
%     result = [0, 0, adj_width, adj_height];

elseif direction == 1
    
    % To get units (x box position):
    screenXpixels = screen(1);
    screenYpixels = screen(2);
    width = image(1);
    height = image(2);
    row_size = rows(1);
    total_images_on_screen = rows(2);
    % Calculate image ratio, and image width:
    img_ratio = height/width;
    adj_width = floor(screenXpixels/row_size);
    adj_height = floor(adj_width * img_ratio);
%     fprintf('\nAdjusted image size: %d, %d \n', adj_width, adj_height);
    adj_height*ceil(total_images_on_screen/row_size);
    img_index = x;
    x_row = rem(img_index-1, row_size) + 1;
    y_row = ceil(img_index / row_size);
%     for i = 1:total_images_on_screen
%         fprintf('y = %i | x = %i\n', ceil(i / row_size), rem(i-1, row_size) + 1)
%     end
    % Get the box for each image:
    result = [ x_row*adj_width-adj_width, y_row*adj_height-adj_height, x_row*adj_width, y_row*adj_height ];
    center = round([ x_row*adj_width-adj_width./2, y_row*adj_height-adj_height./2 ]);
    
elseif direction == 2
    screenXpixels = screen(1);
    screenYpixels = screen(2);
    width = image(1);
    height = image(2);
    row_size = rows(1);
    total_images_on_screen = rows(2);
    % Calculate image ratio, and image width:
    img_ratio = height/width;
    adj_width = floor(screenXpixels/row_size);
    adj_height = floor(adj_width * img_ratio);
%     fprintf('\nAdjusted image size: %d, %d \n', adj_width, adj_height);
    adj_height*(total_images_on_screen/row_size);
    % From x, y coordinates, get img_index:
    x_click = x(1);
    y_click = x(2);
    x_row = ceil(x_click/adj_width);
    y_row = ceil(y_click/adj_height);
    result = (y_row-1)*row_size + x_row;
end
