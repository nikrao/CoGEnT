function data = convertimages(pathname, numsamples, width, height, transmod)
% convertimages - Extracts images patches from all the valid images files in
%                 a directory and adds them to an array that can be used for
%                 learning algorithms.  Uses the even
%
% data = convertimages(pathname, numsamples, width, height)
%
% pathname      - Path with image files to convert
% numsamples    - Number of samples to take from these images
% width         - Width of each sample
% height        - Height of each sample
% transmod      - Transform mod to use (2 for even images)
%
% Returns:
% data          - Matrix with samples (numsamples x (width*height))
%
% JFM   1/9/2001
% Rev:  1/5/2004

data = zeros(numsamples, width*height);

filelist = dir(pathname);
numfiles = 0;

% Find the image files in this directory 
for i = 1:length(filelist);
    % Get filename extension
    [f, ex] = strtok(filelist(i).name, '.');
    
    if(strcmp(ex, '.bmp') == 1 | strcmp(ex, '.jpg') == 1 | ...
       strcmp(ex, '.jpeg') == 1 | strcmp(ex, '.tiff') == 1 | ...
       strcmp(ex, '.tif') == 1 | strcmp(ex, '.pcx') == 1 | ...
       strcmp(ex, '.png') == 1 ) 
   
       % Add transmod # images to data set so that the holdout images 
       % have not even been seen by the dictionary learning algorithms.
   
        [objname, trans_str] = strtok(f, '__');
        trans = sscanf(trans_str(3:length(trans_str)), '%f');
        
        if(mod(trans, transmod) ~= 0)
            numfiles = numfiles + 1;
            imglist(numfiles) = cellstr(filelist(i).name);
        end
        
    end
end

if(numfiles == 0) 
    disp('No images files found (.bmp, .jpg, .jpeg, .tif, .tiff, .pcx, .png)');
    return;
end

n = ceil(numsamples / numfiles);
row = 0;        % Row in data matrix
numzeropatches = 0;

for i = 1:numfiles
    name = strcat(pathname, '\');
    name = strcat(name, char(imglist(i)));
    disp(sprintf('Image: %s', name) );
    img = imread( name );
    [imgheight, imgwidth] = size(img);
    
    for j = 1:n
        row = row + 1;
        if(row > numsamples)
            break;
        end
        
        % Find a random sample of this image
        x = ceil(rand * (imgwidth - width + 1));
        y = ceil(rand * (imgheight - height + 1));
        
        sample = img(y:y+height-1, x:x+width-1);
        
        while( sample == zeros(size(sample)) )
            % If this patchs is all zeros, then try to find another patch
            x = ceil(rand * (imgwidth - width + 1));
            y = ceil(rand * (imgheight - height + 1));
        
            sample = img(y:y+height-1, x:x+width-1);
            
            numzeropatches = numzeropatches + 1;
        end
        
        data(row, :) = reshape(sample', 1, width*height);
        
    end
end

disp(sprintf('Number of all-zero patches found: %d.  These patches were skipped.', numzeropatches));

%save('imgdata', 'data');
