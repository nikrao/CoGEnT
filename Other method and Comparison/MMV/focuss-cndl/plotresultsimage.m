function plotresultsimage(filename, outfilename)
% plotresultsimage   - Loads the results of a trainrd iteration and plots them
%     Meant to be run at the same time are learnrd is running, to split the
%     display and learning processing.  Be careful with variables names in
%     this function, because it loads the workspace from the trainrd.m.
% 
% plotresultsimage(filename, outfilename)
%
%  JFM  1/19/2001
%  Rev: 9/4/2003

true = 1;
false = 0;

printtype = '-dps2';
printext = '.ps';

stop = 0;
plotdegrade = 0;
lastiter = -1;
calc_entropy = true;
img_norm_type = 4;

if(isempty(outfilename) )
   plotsaveresults = 0;
else
   plotsaveresults = 1;
end

while(stop == 0)
 %   [fid, message] = fopen('block', 'r');
 %   if(fid ~= -1)
 %       disp('Could''t open file because of block');
 %       fclose(fid);
 %   else
        %ls
        %keyboard
        %!ls >block
        load(filename);
        workspace = load(params.workspace);
        %!rm block
        
        if(iter ~= lastiter)
            disp(sprintf('Iteration = %d', iter) );
            lastiter = iter;
            
            % ---- Figure 1 ----
            disp('Figure 1...');
            figure(1); clf;
            subplot(4, 1, 1);
            %semilogy(mser); 
            %plot(mser);
            %plot(2:iter, rmse(2:iter));
            semilogy(2:iter, rmse(2:iter));
            
            % If the data was whitened, calculate RMSE in original space
            if(params.whiten == 1)
                % Calculate MSE
                yest = dewhiten*A*x;
                residw = yest - workspace.btrain;
                for k = 1:N
                    errw(k)=norm(residw(:,k))^2;
                end
    
                msew = sum(errw)/(m*N);
                ysigmaw = sqrt(var(reshape(workspace.btrain, m*N, 1)));
                rmsew = sqrt(msew)/ysigmaw;
                fprintf('RMSE/sigma on original (unwhitened) space = %f\n', rmsew);
                
                titlestr = sprintf('Iteration %d  A = (%d x %d), p = %5.3f, RMSE = %f (%f)', ...
                    iter, m, n, p, rmse(iter), rmsew );
            else     
                yest = A*x;
                titlestr = sprintf('Iteration %d  A = (%d x %d), p = %5.3f, RMSE = %f', ...
                    iter, m, n, p, rmse(iter) );
            end
            
            title(titlestr);
            % xlabel('Iteration');
            ylabel('RMSE/\sigma_y');
            setfonts;
            
            subplot(4, 1, 2);
           
            plot(avgdiversity);
            title(sprintf('Average diversity = %f', avgdiversity(iter)) );
            ax = axis;
            ax(4) = n; 
            axis(ax);
          %  xlabel('Iteration');
            ylabel('Average Diversity');
            setfonts;           
 
            
            % Plot p-norm average of x's
            subplot(4, 1, 3);
            
            if(exist('x_pnorm_avg'))
                plot(x_pnorm_avg);
                ylabel('Average ||x||_p');
            end
            
            % Plot A column norm statistics
            subplot(4, 1, 4);
            
            if(exist('Acol_min'))
                hold on;
                plot(Acol_min);
                plot(Acol_max);
                plot(Acol_avg);
                ylabel('Column norm of A');
            end
            
            % -- Format and save figure --
            setfonts;
            set(gcf, 'PaperPosition', [.75 .25 7 10.5]);  % Centered left to right
            
            % Print the date at the bottom right of the page
            h = axes('Position', [0 0 1 1], 'Visible', 'off');
            text(0.8, 0.05, date, 'FontSize', 7);
            
            if(plotsaveresults == 1)
                disp('Printing to file...');
                print(gcf, printtype, sprintf('%s_1%s', outfilename, printext) );
            end
            
            % ---- Figure 4 ----
            % Learned dictionary functions
            disp('Figure 4 (Learned dictionary elements)...');
                        
            figure(4); clf;
            subplot(1,1,1);
            width = params.patchwidth;
            height = params.patchheight;
            outwidth = params.outwidth;
            outheight = params.outheight;
            
            % Sort the columns of A so that the columns resemble each other
            Asort = sortcolumns(A);
            
            if(params.whiten == 1)
                img = sampleimage( (dewhiten * Asort)', width, height, outwidth, outheight, 5);
            else 
                img = sampleimage( Asort', width, height, outwidth, outheight, 5);
            end
            
            image(img);
            colormap(gray(256));
            axis off;
            title(sprintf('Learned dictionary at iteration %d', iter) );
            
            %set(gcf, 'PaperPosition', [.75 .25 7 10.5]);  % Centered left to right
            %set(gcf, 'PaperPosition', [.25 8.0 3 3]);  % Paper position

            setfonts;
           
            % Print the date at the bottom right of the page
           % h = axes('Position', [0 0 1 1], 'Visible', 'off');
           % text(0.8, 0.05, date, 'FontSize', 7);
            
            if(plotsaveresults == 1)
                print(gcf, printtype, sprintf('%s_4%s', outfilename, printext) );
            end

            % ---- Figure 5 ----
            % Original image patches
            disp('Figure 5 (Original image patches)...');
                        
            figure(5); clf;
            subplot(1,1,1);
            width = params.patchwidth;
            height = params.patchheight;
            outwidth = width;
            outheight = height;
            startimage = 768 * 6;
            
            img = sampleimage( (workspace.btrain(:,startimage:startimage+(outwidth*outheight-1)))', ...
                    width, height, outwidth, outheight, img_norm_type);
            
            image(img);
            colormap(gray(256));
            axis off;
            title(sprintf('Original image patches') );
            
            %set(gcf, 'PaperPosition', [.75 .25 7 10.5]);  % Centered left to right
            %set(gcf, 'PaperPosition', [.25 8.0 3 3]);  % Paper position

            setfonts;
           
            % Print the date at the bottom right of the page
           % h = axes('Position', [0 0 1 1], 'Visible', 'off');
           % text(0.8, 0.05, date, 'FontSize', 7);
            
            if(plotsaveresults == 1)
                print(gcf, printtype, sprintf('%s_5%s', outfilename, printext) );
            end
            
            % ---- Figure 6 ----
            % Reconstructed image patches
            disp('Figure 6 (Reconstructed image patches)...');
            
            figure(6); clf;
            subplot(1,1,1);
            width = params.patchwidth;
            height = params.patchheight;
            outwidth = width;
            outheight = height;
            startimage = 768 * 6;
            
            img = sampleimage( (yest(:,startimage:startimage+(outwidth*outheight-1)))', ...
                    width, height, outwidth, outheight, img_norm_type);
            
            image(img);
            colormap(gray(256));
            axis off;
            title(sprintf('Reconstructed image patches') );
            
            %set(gcf, 'PaperPosition', [.75 .25 7 10.5]);  % Centered left to right
            %set(gcf, 'PaperPosition', [.25 8.0 3 3]);  % Paper position

            setfonts;
           
            % Print the date at the bottom right of the page
           % h = axes('Position', [0 0 1 1], 'Visible', 'off');
           % text(0.8, 0.05, date, 'FontSize', 7);
            
            if(plotsaveresults == 1)
                print(gcf, printtype, sprintf('%s_6%s', outfilename, printext) );
            end
            
            
            % ----- Calculate entropy/coding efficiency using eq 25 of Lewicki:1999
            if(calc_entropy == true)
                disp('Entropy...');
            
                nbins = 256;
                
                % Quantize x 
                xq = reshape(x, N*n, 1);
                [ntrash, bin] = hist(xq, nbins);
                binw = bin(2) - bin(1);
                prob = bin + binw/2;
                indx = quantiz(xq, prob);
                bin(nbins+1) = bin(nbins) + binw;
                xq = bin(indx + 1);
                
                xq = reshape(xq, n, N);
                [bits] = entropy(xq, nbins);
                
                mseq = calcmse(A*xq, y);
                rmseq = sqrt(mseq)/ysigma;
	
                
                fprintf('Entropy: # bits >= %f  (bits per pattern)\n', bits);
                fprintf('Entropy: bits/pixel = %f\n', bits / m);
                % fprintf('Entropy2: bits/pixel = %f\n', bits2 * (n/m) );
                fprintf('RMSE = %f  (%f)\n', rmse(iter), rmseq );
                fprintf('Diversity = %f\n', mean(numerosity(x, 1e-4)) );
                
                fout = fopen('entropy.out','a+');
                fprintf(fout, '%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\n', ...
                    params.trialnum, iter, bits/m, rmse(iter), rmseq, avgdiversity(iter), p, nbins, m, n);
                fclose(fout);
            end
            
        end
%    end
    
    % Stop after one interation
    stop = 1;
end
