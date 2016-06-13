function plotresults(filename, outfilename)
% plotresults   - Loads the results of a trainrd iteration and plots them
%     Meant to be run at the same time are learnrd is running, to split the
%     display and learning processing.  Be careful with variables names in
%     this function, because it loads the workspace from the trainrd.m.
% 
% plotresults(filename, outfilename)
%
%  JFM  8/10/2000
%  Rev: 12/2/2004

% B&W for laser printing
%printtype = '-dps2';
%printext = '.ps';

% Encapsulated postscript (EPS) (eps best for Latex)
printtype = '-depsc2';  % Color
%printtype = '-deps2';  % B&W
printext = '.eps';

% Color Windows metafile for poster
%printtype = '-dmeta';
%printext = '.emf';

stop = 0;
lastiter = -1;

if(exist('outfilename') )
   plotsaveresults = 1;
else
   plotsaveresults = 0;
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
        %!rm block
        
        if(iter ~= lastiter)
            disp(sprintf('Iteration = %d', iter) );
            lastiter = iter;
            
            figure(1); clf;
            
            subplot(5, 1, 1);
            %semilogy(mser); 
            %plot(mser);
            plot(2:iter, rmse(2:iter));
            title(sprintf('Iteration %d |A|_F=%5.3f, A = (%d x %d), p = %5.3f, RMSE = %f', iter, norm(A,'fro'), m, n , p, rmse(iter) ));  
           % xlabel('Iteration');
            ylabel('RMSE/\sigma_y');
            setfonts;           
    
            subplot(5, 1, 2);
           
            plot(avgdiversity);
            title(sprintf('Average diversity = %f', avgdiversity(iter)) );
          %  xlabel('Iteration');
            ylabel('Average Diversity');
            setfonts;           
        
            %if(iter > 1)
            %    semilogy(deltax(iter-1, :), '.');
            %end
            %title('Normalized x updates');

        
            subplot(5, 1, 3);
            plot(nummatch);
            ax = axis; ax(4) = n; axis(ax);
           % title(sprintf('Matching columns in A = %d, alpha = %f, gamma = %f, mu = %f', nummatch(iter), alpha, gammaA, mu  ));
            title(sprintf('Matching columns in A = %d, gamma = %f', nummatch(iter), gammaA  ));
            setfonts;           
    
            subplot(5, 1, 4);
            plot(mindist, '.');
            title('Minimum distance');
            axis([1 n 0 1]);
            setfonts;           
    
            subplot(5, 1, 5);
            plot(diversity(iter, :), '.');
            title(sprintf('Diversity (numerosity), %d > %d', sum((diversity(iter,:) > targetdiv)), targetdiv));
            setfonts;           
            
            set(gcf, 'PaperPosition', [.75 .25 7 10.5]);  % Centered left to right
            
            % Print the date at the bottom right of the page
            h = axes('Position', [0 0 1 1], 'Visible', 'off');
            text(0.8, 0.05, date, 'FontSize', 7);
            
            if(plotsaveresults == 1)
                disp('Printing to file...');
                print(gcf, printtype, sprintf('%s_1%s', outfilename, printext) );
            end
            
            % ---- Figure 2 ----
            figure(2); clf;
            
            subplot(5, 1, 1);
            tolerance = params.Atolerance;  %0.01;
            [num, mindist, normf1, matchedframe, index] = compareframe(Aorig, A, tolerance);
            fprintf('Nummatch (within %f) = %f\n', tolerance, num);
            mindist;
            
            % Sparisfy the learned x's to the target diversity
          %  for k = 1:N
          %      xk = x(:,k);
          %      xp(:, k) = (pickhighest(xk, targetdiv))';
          %  end
          
            tolerance = params.Xtolerance;
            [dist, nummatchx1] = comparesolution(xorig, x, index, tolerance);
            fprintf('x solved (tolerance = %f) %d\n', tolerance, nummatchx1);

            plot(dist, '.');
            title(sprintf('Distance ||x_k - x^l_k||, nummatch = %d', nummatchx1));
            setfonts;           
                                  
            
            subplot(5, 1, 2);
            if(exist('nummatchx'))
                plot(nummatchx);
                title(sprintf('Number of x_k matching within %5.2f', params.Xtolerance));        
            end
            setfonts;
            
            
            subplot(5, 1, 3);
         
            semilogy(err, '.');
            title('Residules');
            setfonts;
    
            subplot(5, 1, 4);
            if(iter > 1)
                semilogy(normupdateA);
                titlestr = sprintf('||update A||_F = %f', normupdateA(iter-1) ) ;

                title(sprintf('||update A||_F = %f', normupdateA(iter-1) ) );
            else 
                title(sprintf('||update A||_F ') );
            end
            setfonts;           
            
          %  subplot(5, 1, 4);
          %  plot(lambdaseries);
          %  title(sprintf('Lambda, lambda = %f, lambdamax = %f', lambda, lambdamax) );
          %  setfonts;           
            
            subplot(5, 2, 9);
            plot(x100resid);
            title(sprintf('||Ax100 - b100||  Dist = %f', dist(100) ) );
            setfonts;           

            subplot(5, 2, 10);
            plot(x100norm);
            title(sprintf('||x100|| = %f  Div = %d', x100norm(iter), diversity(iter, 100) ) );
            setfonts;           
    
            set(gcf, 'PaperPosition', [.75 .25 7 10.5]);  % Centered left to right
            
            % Print the date at the bottom right of the page
            h = axes('Position', [0 0 1 1], 'Visible', 'off');
            text(0.8, 0.05, date, 'FontSize', 7);
            
            if(plotsaveresults == 1)
                print(gcf, printtype, sprintf('%s_2%s', outfilename, printext) );
            end
            
            % ---- Figure 3 ----
            figure(3); clf;
            
            subplot(3, 1, 1);

            for i = 1:n
                numzero(i) = sum(x(i,:) <= 1e-3);
            end            
            plot(numzero, '.');
            title(sprintf('Number of times each element of x is zero, avg = %9.1f (expected: %9.1f)', ...
                mean(numzero), N*(1-(avgdiversity(iter)/n))   ) );
            setfonts;
            
            subplot(3, 1, 2);
            
            s100 = sort(abs(x(:,100)));
            plot(s100, '.'); 
            title('Magnitude of elements in x_{100}');
            setfonts;
            
            subplot(3, 1, 3);
            for i = 1:N
                epsilon(i) = sqrt(err(i)) / norm(y(:, i));
                %epsilon(i) = sqrt(err(i)) ;
                %epsilon(i) = norm(y(:,i)) ;

            end
            plot(epsilon, '.');
            title('||y_k - Ax_k|| / ||y_k||');
            
            setfonts;

            set(gcf, 'PaperPosition', [.75 .25 7 10.5]);  % Centered left to right
            
            % Print the date at the bottom right of the page
            h = axes('Position', [0 0 1 1], 'Visible', 'off');
            text(0.8, 0.05, date, 'FontSize', 7);
            
            if(plotsaveresults == 1)
                print(gcf, printtype, sprintf('%s_3%s', outfilename, printext) );
            end            
            
            % ------- Figure 4 (Poster) -----------
            figure(4); clf;
            
            stopiter = 500;
            
            %subplot(3, 1, 1);
            subplot(3, 2, 3);

            plot(avgdiversity);
            ax = axis; ax(2) = stopiter; ax(3) = 0; ax(4) = n; axis(ax);
            title(sprintf('c) Average diversity (n - sparsity) = %5.1f', avgdiversity(iter)) );
          %  xlabel('Iteration');
            ylabel('Non-zero elements of x');
            xlabel('Iteration');
            setfonts;       
            
            %subplot(3, 1, 2);
            subplot(3, 2, 1);

            plot(nummatch);
            ax = axis; ax(2) = stopiter; ax(3) = 0; ax(4) = n; axis(ax);
            title(sprintf('a) Matching columns in A = %d', nummatch(iter) ) );
            ylabel('Number matching');
            %xlabel('Iteration');
            setfonts;           
            
            %subplot(3, 1, 3);
            subplot(3, 2, 2);

            if(exist('nummatchx'))
                plot(nummatchx);
                ax = axis; ax(2) = stopiter; ax(3) = 0; ax(4) = N; axis(ax);
                title(sprintf('b) Matching solutions x_k = %d', nummatchx(iter)) );        
                xlabel('Iteration');
                ylabel('Number matching');
            end            
            setfonts;
            
            % For the poster
            %set(gcf, 'PaperPosition', [.75 .25 7 10.5]);  % Centered left to right
            
            % For the NIPS paper
            set(gcf, 'PaperPosition', [2.5 1 6.0 5]);  % Centered left to right

            if(plotsaveresults == 1)
                print(gcf, printtype, sprintf('%s_4%s', outfilename, printext) );
            end            
            
            % ------- Figure 5 (Poster) -----------
            figure(5); clf;
            
            subplot(3, 1, 1);
            a1 = Aorig(:,5);
            a2 = matchedframe(:,5);
            a1 = a1 / norm(a1);
            a2 = a2 / norm(a2);
            plot(1:m, a1, 'b.', 1:m, a2, 'go');
            ax = axis; ax(2) = m; axis(ax);

            title('Column of A (dot) and matched column of learned A (circle)');
            setfonts;
            
            set(gcf, 'PaperPosition', [.75 .25 7 10.5]);  % Centered left to right

            if(plotsaveresults == 1)
                print(gcf, printtype, sprintf('%s_5%s', outfilename, printext) );
            end            
        
        end
%    end
    
    % Stop after one interation
    stop = 1;
end
