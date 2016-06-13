function Amovie = plotAhistory(Ahist, width, height, outwidth, outheight)
% plotAhistory(Ahistory) - Makes a movie of the development of the A matrix
%             during training.
%
% Ahistory      - Cell array of the learned A matricies 
% width         - Width of each sample
% height        - Height of each sample
% outwidth      - Number of samples across in output img
% outheight     - Number of sample down in output img
%
% JFM   3/6/2001
% Rev:  10/24/2003


numframes = size(Ahist,2);
figure(1);
clf;
axis off;

for i = 1:numframes
    % Plot each frame and add to movie
    
    [img] = sampleimage(Ahist{i}', width, height, outwidth, outheight, 2);
    image(img);
    colormap(gray(256));
    
    title(sprintf('Iteration %d', i));
    
    Amovie(i) = getframe;
end

%movie(Amovie, 1, 10);

