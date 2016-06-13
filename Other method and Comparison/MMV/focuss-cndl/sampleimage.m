function [img] = sampleimage(data, width, height, outwidth, outheight, norm)
% sampleimage -   Converts the rows of matrix data into images of size 
%                 width x height, and arranges them into an image with
%                 outwidth x outheight samples in img.
%
% [img] = sampleimage(data, width, height, outwidth, outheight)
%
%  data     - Data matrix, samples are in each row
%  width    - Width of each sample
%  height   - Height of each sample
%  outwidth - Number of samples across in output img
%  outheight- Number of sample down in output img
%  norm     - = 1 to normalize (for input real-valued images, 0.0 -> 128)
%             = 2 to remove mean (for A with prewhitening)
%             = 3 no normalization, multiply by 256, add 128
%             = 4 no normalization, multiply by 256
%             = 5 to center at 0, and normalize max of whole data
% Returns:
%  img      - Composite image (matrix)
%
% JFM   1/10/2001
% Rev:  9/4/2003

index = 1;
[datarows, datacols] = size(data);
y = 1;

img = zeros(outheight * height + (outheight-1)*1, outwidth * width + (outwidth-1) * 1);
img = img + 255;

max_data = max(max(abs(data)));

for row = 1:outheight
    x = 1;
    for col = 1:outwidth
        if(index > datarows)
            break;
        end
                
        sample = real(reshape(data(index, :), width, height))';
        
        if(norm == 1)
            % Normalize the sample to use the full range of grays
            mx = max(max(abs(sample)));
            if(mx == 0) 
                mx = 1;
            end
                            
            sample = 127*(sample ./ mx) + 128;
            
        elseif(norm == 2)
            sample = sample - mean(mean(sample));
            
            mx = max(max(abs(sample)));
            if(mx == 0) 
                mx = 1;
            end
                            
            sample = 127*(sample ./ mx) + 128;
            
        elseif(norm == 3)
            sample = 256 * sample + 128;
            
        elseif(norm == 4)
            sample = 256 * sample;            
         
        elseif(norm == 5)
            sample = 127 *(sample ./ max_data) + 128;
            
        end
        
        img(y:y+height-1, x:x+width-1) = sample;
       
        index = index + 1;
        x = x + width + 1;
    end
    
    y = y + height + 1;
end

img = uint8(img);