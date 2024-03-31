function ellipsePixels = createEllipse(imageSizeX, imageSizeY)
% Create a logical image of an ellipse with specified
% semi-major and semi-minor axes, center, and image size.
% First create the image.

imageSizeX = imageSizeX;
imageSizeY = imageSizeY;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

% Next create the ellipse in the image.
centerX = imageSizeX/2 + 1;
centerY = imageSizeY/2 + 1;
radiusX = imageSizeX/2 + 1;
radiusY = imageSizeY/2 + 1;
ellipsePixels = (rowsInImage - centerY).^2 ./ radiusY^2 ...
    + (columnsInImage - centerX).^2 ./ radiusX^2 <= 1;
	
% ellipsePixels is a 2D "logical" array.
% Now, display it.
% image(ellipsePixels) ;
% colormap([0 0 0; 1 1 1]);
% title('Binary image of a ellipse', 'FontSize', 20);
end
