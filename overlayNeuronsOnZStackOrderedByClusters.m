function out = overlayNeuronsOnZStackOrderedByClusters(coors,allRawTracesNormalized,validNeuronIDs,imageName,clusters,uniqueColors,out);
    %this code will loop through every image in the stack and generate a
    %stacked line plot next to the movie in order to visualise activity over
    %time side-by-side with the tif.
    %this will colour each line by the neuron's assigned clusters, useful
    %to validate if all neurons in a cluster are behaving the same way
    Coor=coors;
    timeCourseLength=size(out.C,2);

    
    t = Tiff(imageName,'r');
    InfoImage=imfinfo(imageName);
    offsets=length(InfoImage);
    
    largestCluster=0;
    for clusterIndex=1:length(clusters)
        neuronsInThisCluster=length(clusters{clusterIndex});
        if neuronsInThisCluster>largestCluster
            largestCluster=neuronsInThisCluster;
        end
    end
    cellCrops=zeros(50*length(clusters),50*largestCluster,3,'uint8');
    
    yScalar=0;
    xScalar=0;
    timestamp=datestr(datetime,'yyyymmddTHHMMSS');
    imageName=strcat('outputMovieClustered',timestamp,'.tif');
    %adjust to lower temporal resolution if image size too large or small
    %for i= 1:offsets
    for i= 1:50:offsets
    %for i= 1:25

        %subimage=imread(t,'Index',i,'Info',InfoImage);
        t.setDirectory(i);
        subimage=t.read();
        tempFig=figure();
        thisFrameImage=imagesc(subimage);
        hold on;
        
        if 1==i
            thisFrame = getframe(gca);
            originalSize=size(subimage);
            originalyHeight=originalSize(1);
            originalxHeight=originalSize(2);
            framesize=size(thisFrame.cdata);
            frameyHeight=framesize(1);
            framexHeight=framesize(2);
            yScalar=frameyHeight/originalyHeight;
            xScalar=framexHeight/originalxHeight;
        end
        
        
        for clusterIndex=1:length(clusters)
            neuronsInThisCluster=clusters{clusterIndex};
            clusterCenters=zeros(2,length(neuronsInThisCluster));
            activationsInCluster=0;
            thisClusterInterpeakTimes=[];
            currentOffset=0;
            for j=1:length(neuronsInThisCluster)
                thisNeuron=neuronsInThisCluster(j);
                plot(Coor{thisNeuron}(1,:),Coor{thisNeuron}(2,:), 'color',uniqueColors(clusterIndex,:),'lineWidth',1);
                minX=min(Coor{thisNeuron}(1,:));
                minY=min(Coor{thisNeuron}(2,:));
                maxX=max(Coor{thisNeuron}(1,:));
                maxY=max(Coor{thisNeuron}(2,:));

                minX=minX*xScalar;
                minY=minY*yScalar;
                
                thisFrame = getframe(gca);
                thisFrame=thisFrame.cdata;
                thisCellCrop=imcrop(thisFrame,[max(minX-5,0),max(minY-5,0),50,50]);
                cellCrops((((clusterIndex-1)*50)+1):((((clusterIndex-1)*50)+1)+49),(((j-1)*50)+1):((((j-1)*50)+1)+49),:)=imresize(thisCellCrop,[50,50]);
                
                
                clusterColors(:,clusterIndex)=uniqueColors(clusterIndex,:)';
                clusterCenters(:,j)=[Coor{thisNeuron}(1,1);Coor{thisNeuron}(2,1)];
                thisNeuronTimeCourse=out.S(thisNeuron,:);
                spikesInThisTimecourse=thisNeuronTimeCourse>0;
                totalSpikesInThisTimescourse=sum(spikesInThisTimecourse);
                activationsInCluster=activationsInCluster+totalSpikesInThisTimescourse;
                lastPeakTime=0;
                thisNeuronInterpeakTimes=[];
                for k=1:length(spikesInThisTimecourse)
                    if spikesInThisTimecourse(k)==1
                        if(lastPeakTime==0)
                            lastPeakTime=k;
                        else
                           lastPeakTime=k;
                        end
                    end
                end
            end
        end
        close(tempFig);
        cellCropHeight=size(cellCrops,1);
        set(gcf,'color','k');
        F1 = getframe(gcf);
        %movieData=F1.cdata;
        
        
        
        %figure('visible','off');
        tempFig=figure();
        
        
        hold on;
        neuronIndex=1;
        for clusterIndex=1:length(clusters)
            neuronsInThisCluster=clusters{clusterIndex};
            for thisNeuronIndex=1:length(neuronsInThisCluster)
                plot3([1:max(2,i)],repmat(neuronIndex,max(2,i)),allRawTracesNormalized(1:max(2,i),neuronsInThisCluster(thisNeuronIndex)),'lineWidth',2,'Color',uniqueColors(clusterIndex,:));
                neuronIndex=neuronIndex+1;
                if thisNeuronIndex==length(neuronsInThisCluster)
                    neuronIndex=neuronIndex+3;
                end
            end
        end
        view(0,88);
        
        set(gcf,'color','k');
        set(gca,'color','k');
        xlim([1,offsets]);

        F2 = getframe(gcf);
        graphData = F2.cdata;
        close(tempFig);

        maxXDimension=size(subimage,2)+size(graphData,2);
        maxYDimension=max(size(subimage,2),size(graphData,2));
        
        %arrange the movie data and the plot data side by side in the
        %correct orientation
       	movieData=imresize(cellCrops,size(graphData,1)/cellCropHeight);
        imdata=[movieData,graphData];
        imdata=imresize(imdata,[250,NaN]);
        if i==1
            options.color = true;
            options.compress='lzw';
            saveastiff(imdata, imageName,options);
        else
            clear options;
            options.compress='lzw';
            options.color = true;            
            options.append = true; 
            timestamp=datestr(datetime,'yyyymmddTHHMMSS');
            saveastiff(imdata, imageName,options);
        end
    end

end
