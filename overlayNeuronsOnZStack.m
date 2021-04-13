function out = overlayNeuronsOnZStack(coors,allRawTracesNormalized,validNeuronIDs,imageName);
    %this code will loop through every image in the stack and generate a
    %stacked line plot next to the movie in order to visualise activity over
    %time side-by-side with the tif.
    
    %requires stackedPlot from file exchange
    % www.mathworks.com/matlabcentral/fileexchange/24368-stacked-plot/
    
    t = Tiff(imageName,'r');
    InfoImage=imfinfo(imageName);
    offsets=length(InfoImage);
    
    %loop through every image in the stack
    for i= 1:offsets
        
        %subimage=imread(t,'Index',i,'Info',InfoImage);
        t.setDirectory(i);
        subimage=t.read();
        figure('visible','off');
        thisFrameImage=imagesc(subimage);
        hold on;
        for k=1:length(coors)
            plot(coors{k}(1,:),coors{k}(2,:), '-w','lineWidth',1);
            text(max(coors{k}(1,:))-(max(coors{k}(1,:))-min(coors{k}(1,:))),min(coors{k}(2,:))-10,num2str(k),'Color','w','FontSize',12);
        end
        set(gcf,'color','k');
        F1 = getframe(gcf);
        movieData=F1.cdata;
        
        figure('visible','off');
        h=stackedplot(allRawTracesNormalized(1:max(2,i),:),3,1,2,'faceColor','none','lineWidth',2);
        set(gca,'GridColor','w');
        set(gca, 'GridAlpha',0.25);
        set(gca, 'XColor','w');
        xticks([1:2:length(coors)]);
        
        set(gcf,'color','k');
        set(gca,'color','k');
        ylim([1,offsets]);
        grid on;
        colormap('parula');
        ylabel('time');
        view(-90.01,85);
        %set(gcf,'color','w');
        F2 = getframe(gcf);
        %[X, ~] = frame2im(F);        
        graphData = F2.cdata;
       
        maxXDimension=size(subimage,2)+size(graphData,2);
        maxYDimension=max(size(subimage,2),size(graphData,2));
        
        
        imdata=[movieData,graphData];
        if i==1
            options.color = true;
            options.compress='lzw';
            saveastiff(imdata, 'outputMovie.tif',options);
        else
            clear options;
            options.compress='lzw';
            options.color = true;            
            options.append = true;            
            saveastiff(imdata, 'outputMovie.tif',options);
        end
    end

end
