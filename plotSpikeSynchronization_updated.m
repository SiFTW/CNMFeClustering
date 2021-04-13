%synchronized activity figure
%fileName='neuronOut (250frm) M+C isoCtrl 05-10#2.mat';
%fileName='neuronOut M+C RTT (sync) 07-06#3.mat';
fileName='neuronOut M+C RTT (sync) 07-06#3';
load(fileName);

%syncThreshold=0.4;
syncThreshold='auto';
numberOfShufflesToGenerate=1000;
percentileForShuffling=95;

mySpikes= out.S;
mySpikes(mySpikes>0)=1;


synchronizationShuffled=zeros(1,size(mySpikes,2),numberOfShufflesToGenerate);

if strcmp(syncThreshold,'auto')
    %perform spike shuffling to determine the threshold for xcorr activity.
    shuffledTraces=zeros(size(mySpikes,1),size(mySpikes,2),numberOfShufflesToGenerate);
    for shuffle=1:numberOfShufflesToGenerate
        for i=1:size(shuffledTraces,1)
            shuffledTraces(:,i,shuffle)=mySpikes(randperm(size(mySpikes,1)),i);
        end
    end
    %now calculate the synchronization of the shuffled spikes over time
    for j=i:numberOfShufflesToGenerate
        for i=1:size(synchronizationShuffled,2)
            synchronizationShuffled(1,i,j)=sum(sum(shuffledTraces(:,max(1,i-2):min(i+2,size(shuffledTraces,2)),j)));
        end
    end
    synchronizationShuffled=synchronizationShuffled./size(mySpikes,1);
    synchronizationShuffled(synchronizationShuffled>1)=1;
    synchronizationThresh=prctile(synchronizationShuffled(:),percentileForShuffling);
else
    synchronizationThresh=syncThreshold;
end

%R=abs(R)
figure();
surf(mySpikes);
view(2);
set(gcf,'color','w');
saveas(gcf,'spikes.fig')


figure;
subplot(6,1,[1,2]);
mySpikes= out.S;
mySpikes(mySpikes>0)=1;
h=surf(mySpikes);
set(h,'linestyle','none');
view(2);
set(gcf,'color','w');
set(gca,'yTick',[]);
set(gca,'xTick',[]);


subplot(6,1,3);
synchronization=zeros(1,size(mySpikes,2));
for i=1:length(synchronization)
    synchronization(i)=sum(sum(mySpikes(:,max(1,i-2):min(i+2,size(mySpikes,2)))));
end
synchronization=synchronization./size(mySpikes,1);
synchronization(synchronization>1)=1;
plot(synchronization,'k','lineWidth',2);
hold on;
plot([0,length(synchronization)],[synchronizationThresh,synchronizationThresh],'r','lineWidth',1);
set(gca,'yTick',[0,1]);
ylim([0,1]);
set(gca,'xTick',[]);
framesSynced=sum(synchronization>synchronizationThresh);
peaksAbove=bwconncomp(synchronization>synchronizationThresh);
SyncEvents=peaksAbove.NumObjects;
xlabel(sprintf('frames synced:%d, num of sync events: %d',framesSynced,SyncEvents));

subplot(6,1,4);
h=surf([synchronization;synchronization]);
set(h,'linestyle','none');
view(2);
set(gca,'yTick',[]);
set(gca,'xTick',[]);
caxis([0,1]);

subplot(6,1,[5,6]);
stackedplot(out.C',3,1);
grid off;
colormap('parula');
ylabel('time');
view(-90,75);
set(gca,'yTick',[1:50:251]);
set(gca,'yTickLabel',[0:50:250]);
set(gca,'xTick',[]);
set(gca,'zTick',[]);
