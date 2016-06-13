figure; hold on
for i=1:5
    for j=1:5
        for k=1:5
            z=rand(1,1);
            if z<=.5
                scatter3(i,j,k,'k o')
            else
                scatter3(i,j,k,'r .')
            end
        end
    end
end