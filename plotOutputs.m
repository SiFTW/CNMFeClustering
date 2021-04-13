function out =  plotOutputs()
%used to generate figures in Samarasinghe et al.

%this is the output file name
filename='output';

%the type of correlation to perform
correlationType='xcorr'; %xcorr,corr,FCA

%the file name of the tif stack used to make the movies
imageName='images.tif';

%the output file from CNMFE
%https://github.com/zhoupc/CNMF_E
load('neuronOut M+C RTT (sync) 07-06#3.mat');

%theshold linkage parameter, at which the dendogram is split into clusters.
linkageParam=1.2;

%how many times should spikes be shuffled to get the maximum correlation of
%shufled neurons
numberOfShufflesToGenerate=1000;

%If using xcorr, and a neuron only spikes once, the xcorr will be perfect
%so lets drop there neurons
if(strcmp(correlationType,'xcorr'))
    validNeuronIds=find(out.P.neuron_sn>0);
    onespikeNeurons=[];
    for i=1:length(validNeuronIds)
        if length(find(out.S(validNeuronIds(i),:)))<=1
            onespikeNeurons=[onespikeNeurons,validNeuronIds(i)];
        end
    end
    validNeuronIds(onespikeNeurons)=[];
else
    validNeuronIds=find(out.P.neuron_sn>0);
end


%keep track of neuron IDs
neuronIDs=[1:length(validNeuronIds)];
numberOfNeurons=length(validNeuronIds);
%the length of the timecourse is determined from CNMFE output
timeCourseLength=size(out.C,2);
%create some variables for storing neuron traces
%just the traces straight from CNMFE
allRawTraces=zeros(numberOfNeurons,timeCourseLength);
%the traces normalised to the max of this neurons trace
allRawTracesNormalized=zeros(numberOfNeurons,timeCourseLength);

%we need interpeak times to calculate burstiness and memory, initiate some
%empty arays
meanInterPeakTimePerNeuron=zeros(numberOfNeurons,1);
allRawInterpeakTimes={};
allInterpeakTimesWithoutClusters={};
arrayOfBustines=zeros(numberOfNeurons,1);
arrayOfMemoryWithCluster=zeros(numberOfNeurons,1);
arrayOfMemoryWithoutCluster=zeros(numberOfNeurons,1);

%go through each neuron
for i=1:length(validNeuronIds)
    thisID=validNeuronIds(i);
    allRawTraces(i,:)=out.C(thisID,:);
    %divide by maximum value in trace to get a normalised timecourse
    allRawTracesNormalized(i,:)=allRawTraces(i,:)./max(allRawTraces(i,:));
    %get the spiking data from CNMFE output
    spikeCourse=out.S(thisID,:);
    %initialise an empty array of interpeak times
    thisNeuronPeakTimes=[];
    lastPeak=0;
    %go through every spike
    for k=1:length(spikeCourse)
        %if there is a spike
        if(spikeCourse(k)>0)
            %if not the first spike
            if(lastPeak~=0)
                %append how long from the previous spike to this one
                thisNeuronPeakTimes=[thisNeuronPeakTimes,k-lastPeak];
            end
            %the last spike occured at this time
            lastPeak=k;
        end
    end
    
    %as long as there are enough interpeak times to calculate memory and
    %burstiness for each cluster of neurons lets do it
    if(length(thisNeuronPeakTimes)>2)
        allRawInterpeakTimes{i}=thisNeuronPeakTimes;
        allInterpeakTimesWithoutClusters{i}=thisNeuronPeakTimes(thisNeuronPeakTimes>1);
        meanInterPeakTimeThisNeuron=mean(thisNeuronPeakTimes);
        meanInterPeakTimePerNeuron(i)=meanInterPeakTimeThisNeuron;
        
        stdInterPeakTimeThisNeuron=std(thisNeuronPeakTimes);
        arrayOfBustines(i)=(stdInterPeakTimeThisNeuron-meanInterPeakTimeThisNeuron)/...
            (stdInterPeakTimeThisNeuron+meanInterPeakTimeThisNeuron);
        [arrayOfMemoryWithCluster(i),~]= corr(thisNeuronPeakTimes(1:end-1)',thisNeuronPeakTimes(2:end)','type', 'Spearman');
        if(length(allInterpeakTimesWithoutClusters{i})>2)
            [arrayOfMemoryWithoutCluster(i),~]= corr(allInterpeakTimesWithoutClusters{i}(1:end-1)',allInterpeakTimesWithoutClusters{i}(2:end)','type', 'Spearman');
        end

    else
        % not enough data for calculations
        arrayOfBustines(i)=NaN;
        arrayOfMemoryWithCluster(i)=NaN;
        arrayOfMemoryWithoutCluster(i)=NaN;
    end

end

%%% uncomment to make zstack image
% requires stacked plot from  MATLAB File Exchange
% www.mathworks.com/matlabcentral/fileexchange/24368-stacked-plot/

% overlayNeuronsOnZStack(out.Coor,allRawTracesNormalized',validNeuronIds,imageName);

%%%%%% calculate the overall interpeak times for all neurons, rather than
%%%%%% for clusters above
%loop through the time course
%caculate time since last spike
populationSpikes=zeros(size(out.S,2));
populationInterpeakTimes=[];
lastSpikeTime=0;
for i=1:size(out.S,2)
    spikeAtTimepoint=0;
    for j=1:size(out.S,1)
        if(out.S(j,i)>0)
            spikeAtTimepoint=1;
        end
    end
    populationSpikes(i)=spikeAtTimepoint;
    if(spikeAtTimepoint==1)
        timeSinceLastSpike=i-lastSpikeTime;
        populationInterpeakTimes(end+1)=timeSinceLastSpike;
        lastSpikeTime=i;
    end
end
%calculate burstiness and memory for the entire population
meanInterPeakTimePopulation=mean(populationInterpeakTimes);
stdInterPeakTimePopulation=std(populationInterpeakTimes);
burstinessPopulation=(stdInterPeakTimePopulation-meanInterPeakTimePopulation)/...
(stdInterPeakTimePopulation+meanInterPeakTimePopulation);
memoryPopulation= corr(populationInterpeakTimes(1:end-1)',populationInterpeakTimes(2:end)','type', 'Spearman');


%display bustiness as a histogram
figure();
histogram(arrayOfBustines);
xlabel('Burstiness Per Neuron');
ylabel('count')

%display memory as a histogram
figure();
histogram(arrayOfMemoryWithCluster);
xlabel('Memory Per Neuron');
ylabel('count')


%Save a few figures showing the traces
figure();
% requires stacked plot from  MATLAB File Exchange
% www.mathworks.com/matlabcentral/fileexchange/24368-stacked-plot/
stackedplot(allRawTracesNormalized',3,1);
grid off;
colormap('parula');
ylabel('time');
view(-90,85);
saveas(gcf,'traces.fig')
figure();
stackedplot(allRawTracesNormalized',5,1);
xlabel('time');
saveas(gcf,'topviewtraces.fig')

%histograms of interpeak times.
figure();
histogram(meanInterPeakTimePerNeuron((meanInterPeakTimePerNeuron>0)));
xlabel('mean interpeak times per neuron (s)');
ylabel('count')
saveas(gcf,'interpeak.fig')

% plot the neurons on the image, with each neurons boundary plotted.
ax1=figure();
imagesc(out.Cn, [0, 1]);
colormap(ax1,'gray');
hold on;
axis off;
for k=1:length(validNeuronIds)
    plot(out.Coor{validNeuronIds(k)}(1,:),out.Coor{validNeuronIds(k)}(2,:), '-b','lineWidth',1);
end
saveas(gcf,'selected.fig')


%perform spike shuffling to determine the threshold for xcorr activity.
shuffledTraces=zeros(size(allRawTracesNormalized,1),size(allRawTracesNormalized,2),numberOfShufflesToGenerate);
for shuffle=1:numberOfShufflesToGenerate
    for i=1:size(shuffledTraces,1)
        shuffledTraces(:,i,shuffle)=allRawTracesNormalized(randperm(size(allRawTracesNormalized,1)),i);
    end
end
R=zeros(size(allRawTracesNormalized,1));
Rshuffles=zeros(size(R,1),size(R,2),numberOfShufflesToGenerate);
if(strcmp(correlationType,'corr'))
    R=abs(corrcoef(allRawTracesNormalized'));
elseif(strcmp(correlationType,'xcorr'))
    for i=1:size(R,1)
        for j=1:size(R,1)
            R(i,j)=max(xcorr(allRawTracesNormalized(:,i),allRawTracesNormalized(:,j),'coeff'));
            for k=1:numberOfShufflesToGenerate
                Rshuffles(i,j,k)=max(xcorr(shuffledTraces(:,i,k),shuffledTraces(:,j,k),'coeff'));
            end
        end
    end
    %R=abs(xcorr(allRawTracesNormalized','coeff'));
elseif(strcmp(correlationType,'fca'))
    %FCA code is untested
    fileID=fopen('FCAInputFile.txt','w');
    for i=1:size(allRawTracesNormalized,1)
        thisTimeCourse=out.S(i,:);
        spikeTimes=find(thisTimeCourse>0);
        for j=1:length(spikeTimes)
            fprintf(fileID,'%d %d\n',i,spikeTimes(j));
        end

    end
    fclose(fileID);
end


%lets calculate significant pairs vs the shuffled correlations.
sigVals=zeros(size(R));
for i=1:size(R,1)
    for j=1:size(R,1)
        sigVals(i,j)=ttest2(Rshuffles(i,j,:),R(i,j));
    end
end
numberOfPairs=sum(sum(triu(sigVals,1)));
denominator=sum(sum(triu(ones(size(sigVals)),1)));
csvwrite('numberOfPairs.csv',numberOfPairs);
csvwrite('percentageOfPairs.csv',numberOfPairs./denominator);
upperTriR1=triu(R,1);
csvwrite('averageCorrCoeff.csv',mean(upperTriR1(:)));

if(strcmp(correlationType,'fca'))
    %additional untested FCA code
    clustering=dlmread('FCAOutputclust.txt',',');
    distances=dlmread('FCAOutputdist.txt',',');
    tree=[clustering,distances];
    tree(:,3)=cumsum(tree(:,3));
    figure;
    colormap(gcf,jet) 
    [H,T,Outperm]=dendrogram(tree,0,'ColorThreshold',13);  
else
    %R=abs(R)
    value=sum(R(:)>0.3);
    value = value - size(R,1);
    totalSquares=(length(R(:))-size(R,1))/2;
    value=value/2;
    percentPositive=(value/totalSquares)*100
    figure();
    mySpikes= out.S;
    mySpikes(mySpikes>0)=1;
    surf(mySpikes);
    view(2);
    set(gcf,'color','w');
    saveas(gcf,'spikes.fig')

    %Percent of >30% correlated spikes
    %if we can't do FCA we'll use an abitrary cutoff for R2 but this is
    %output and not used for rest of analysis
    R2=abs(corrcoef(mySpikes'));
    value2=sum(R2(:)>0.3);
    value2 = value2 - size(R2,1);
    totalSquares=(length(R2(:))-size(R2,1))/2;
    value2=value2/2;
    percentPositive3=(value2/totalSquares)*100

    %cluster the neurons using their correlations
    cgo= clustergram(R,'colormap','jet','ShowDendrogram',true,'symmetric',false);
    %use the previous linkage parameter
    set(cgo,'Linkage','complete','Dendrogram',linkageParam)
    set(gcf,'color','w');
    caxis([0,1]);

    figure;
    tree=linkage(R,'complete');
    figure;
    colormap(gcf,jet)
    %output just the dendrogram
    [H,T,Outperm]=dendrogram(tree,0,'Reorder',optimalleaforder(tree,pdist(R)),'ColorThreshold',linkageParam);
end
%save the dendrogram
set(H,'LineWidth',3)
colorArray=zeros(3,length(H));
set(gcf,'color','w');
daspect([2 1 1]);
box on;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
saveas(gcf,'dendrogram.png');

%loop through the lines in the dendrogram and get the colours, use this
%colour info to assign each neuron to a cluster
for i=1:length(H)
    thisLine=H(i);
    %thisColor=thisLine.Color;
    thisColor=get(H(i),'Color');
    %thisYData=thisLine.YData;
    thisYData=get(H(i),'YData');
    %thisXData=thisLine.XData;
    thisXData=get(H(i),'XData');
    for j=1:length(thisYData)
        %if this is a leaf node
        if thisYData(j)==0
            thisXIndex=thisXData(j);
            leafIndex=Outperm(thisXIndex);
            %here we store the colour of each neuron
            colorArray(:,leafIndex)=thisColor;
        end
    end
end
%uniqueColours will tell us how many clusters we have and what colour each
%one is
[uniqueColors,~,colorIndex]=unique(colorArray','rows');
clusters=cell(1,length(uniqueColors));
for i=1:length(colorIndex)
    thisCluster=colorIndex(i);
    if isempty(clusters{thisCluster})
        clusters{thisCluster}=i;
    else
        clusters{thisCluster}(end+1)=i;
    end
end

%let's plot the clusters spatially back on the tiff image
%this time we'll plot them in colour with the neuron borders in each colour
%and the image in grayscale.
figure;
t = Tiff(imageName,'r');
t=t.read();
thisImage=imagesc(t);
hold on;
colormap('gray');

%we're going to store some metrics for each cluster
%initialise those variables here
spatialExtentOfClusters=zeros(1,length(clusters));
totalSpatialExtentOfClusters=zeros(1,length(clusters));
neuronsInCluster=zeros(1,length(clusters));
clusterColors=zeros(3,length(clusters));
averagePairwiseDistances=zeros(1,length(clusters));
maximumPairwiseDistances=zeros(1,length(clusters));
averageActivationsPerFrame=zeros(1,length(clusters));
timeCourseLength=length(out.S(1,:));
totalActivationsInCluster=zeros(1,length(clusters));
averageInterpeakTimes=zeros(1,length(clusters));
nonZeroElements=[];

% we only care about non-empty clusters
% drop the rest
for i=1:length(clusters)
    if(~isempty(clusters{i})&&length(clusters{i})>2)
        nonZeroElements=[nonZeroElements,i];
    end
end
clusters=clusters(nonZeroElements);

%go through each cluster and plot its neurons on the image in the right color
for i=1:length(clusters)
    neuronsInThisCluster=clusters{i};
    clusterCenters=zeros(2,length(neuronsInThisCluster));
    activationsInCluster=0;
    thisClusterInterpeakTimes=[];
    %for every neuron in the cluster
    for j=1:length(neuronsInThisCluster)
        thisNeuron=neuronsInThisCluster(j);
        %plot its outline in the right colour
        plot(Coor{thisNeuron}(1,:),Coor{thisNeuron}(2,:), 'color',uniqueColors(i,:),'lineWidth',1);
        clusterColors(:,i)=uniqueColors(i,:)';
        clusterCenters(:,j)=[Coor{thisNeuron}(1,1);Coor{thisNeuron}(2,1)];
        thisNeuronTimeCourse=out.S(thisNeuron,:);
        spikesInThisTimecourse=thisNeuronTimeCourse>0;
        totalSpikesInThisTimescourse=sum(spikesInThisTimecourse);
        activationsInCluster=activationsInCluster+totalSpikesInThisTimescourse;
        lastPeakTime=0;
        thisNeuronInterpeakTimes=[];
        %lets caluclate interpeak times for this cluster
        for k=1:length(spikesInThisTimecourse)
            if spikesInThisTimecourse(k)==1
                if(lastPeakTime==0)
                    lastPeakTime=k;
                else
                   thisNeuronInterpeakTimes(end+1)=k-lastPeakTime;
                   lastPeakTime=k;
                end
            end
        end
        thisClusterInterpeakTimes=[thisClusterInterpeakTimes,thisNeuronInterpeakTimes];
    end
    xPoints=clusterCenters(1,:);
    yPoints=clusterCenters(2,:);
    
    %get the boundary of the cluster
    [boundaryPoints,boundaryArea]=boundary(xPoints',yPoints');
    totalSpatialExtentOfClusters(i)=boundaryArea;
    
    neuronsInCluster(i)=length(neuronsInThisCluster);
    
    %how big is the cluster normalsied by the number of neurons in it
    spatialExtentOfCluster=boundaryArea/length(neuronsInThisCluster);
    spatialExtentOfClusters(i)=spatialExtentOfCluster;
    
    %plot the cluster area as a patch on the image
    patch('XData',xPoints(boundaryPoints),'YData',yPoints(boundaryPoints),'FaceColor',uniqueColors(i,:),'FaceAlpha',.25,'EdgeColor',uniqueColors(i,:),'LineWidth',2);
    
    %mean pairwise distance between neurons in the cluster
    averagePairwiseDistance=mean(pdist([xPoints',yPoints']));
    averagePairwiseDistances(i)=averagePairwiseDistance;
    %max pairwise distance between neurons in the cluster
    maximumPairwiseDistance=max((pdist([xPoints',yPoints'])));
    maximumPairwiseDistances(i)=maximumPairwiseDistance;
    %spikes in the cluster
    totalActivationsInCluster(i)=activationsInCluster;
    %spikes per frame in the cluster
    averageActivationsPerFrame(i)=activationsInCluster/(timeCourseLength*length(neuronsInThisCluster));
    %mean interpeak time in the cluster
    averageInterpeakTimesInCluster=mean(thisClusterInterpeakTimes);
    averageInterpeakTimes(i)=averageInterpeakTimesInCluster;
    
end
%save this image of all the cluster areas colour coded
saveas(gcf,'clusterLocations.png');


%let's plot all neuron traces but coloured and ordered by clustering
%recomment this to plot a video of the neurons, ordered by clusters, side
%by side with their clustered traces.
%overlayNeuronsOnZStackOrderedByClusters(out.Coor,allRawTracesNormalized',validNeuronIds,imageName,clusters,uniqueColors,out);


%Below code counts the number of neurons "S" that are spiking at any given
%time.  PercentSpikeCorr is the percentage of the total that are spiking.
M=mean(mySpikes,2);
Smean=sum(M);
Si=size(mySpikes,1)
Avg=((Smean/Si)/0.3)
S=sum(mySpikes);
PercentSpikeCorr=S/Si;

%This is to calculate burtsiness by the Schleiss and Smith methodology, M
%is defined above
StDev=std(mySpikes,0,2);
A1=StDev-M;
B1=StDev+M;
Burstiness=A1./B1;


csvwrite('myfileBursty.csv',Burstiness);
csvwrite('myfileR.csv',R);
csvwrite('myfileR2.csv',R2);
%csvwrite('myfileIntP.csv', interpeakTime);
csvwrite('myfileNormtrace.csv',allRawTracesNormalized);
csvwrite('myfileSpikes.csv',mySpikes);

csvwrite('meatInterpeakTimePerNeuron.csv',meanInterPeakTimePerNeuron);
csvwrite('burstiness.csv',arrayOfBustines);
csvwrite('memory.csv',arrayOfMemoryWithCluster);
csvwrite('burstinessAndMemoryPopulation.csv',[burstinessPopulation,memoryPopulation]);
csvwrite('averagePairwiseDistances.csv',averagePairwiseDistances);
csvwrite('spatialExtentOfClusters.csv',spatialExtentOfClusters);
csvwrite('maximumPairwiseDistances.csv',maximumPairwiseDistances);
csvwrite('averageActivationsPerFrame.csv',averageActivationsPerFrame);
csvwrite('clusterColors.csv',clusterColors);
csvwrite('totalSpatialExtentOfClusters.csv',totalSpatialExtentOfClusters);
csvwrite('neuronsInCluster.csv',neuronsInCluster);
csvwrite('totalActivationsInCluster.csv',totalActivationsInCluster);
csvwrite('averageInterpeakTimes.csv',averageInterpeakTimes);

end
