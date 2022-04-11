%% Clear previous data
clc; close; clear
clear all
filename=' Trackit 83 SRF_Probe_16_2_20.mat'
load (filename(1:end))
AnalyseMin=1;
timesec=AnalyseMin*60;
timesec=59
realworldvalue=92;
fps=30;
filt=1;
if filt>1
filteredcom=movmedian(filteredcom,filt);
filterednose=movmedian(filterednose,filt);
end
I=(BackgroundCropped);
se = strel('disk',5);
background = imdilate(I,se);
maze2 = edge(background,'Canny');
maze2 = bwareaopen(maze2,600);
maze=(imfill(maze2,'holes'));
maze = imerode(maze,se);
I(~maze) = 255;
Cmaze = regionprops(maze,'Centroid','MajorAxisLength','Area');
conversionfactor=realworldvalue/Cmaze.MajorAxisLength;
% Distance Travelled
r =timesec*fps;% for reshaping by frame rate 
dstnc=[];
xt=filterednose(:,1);
yt=filterednose(:,2);
xt=filterednose(:,1);
yt=filterednose(:,2);
f=length(filterednose);
for i=1:f-1
d = sqrt((xt(i+1)-xt(i))^2+(yt(i+1)-yt(i))^2); % calculate distance
dstnc=[dstnc;d];
end
Frames=dstnc(1:r);
DT = (sum(reshape(Frames,fps,[]),1))*conversionfactor;
% escape hole
imshow(BackgroundCropped)
L = false(size(I));
hold on
FirstEscape = drawpoint('Color','k');
mask=false(size(I));
coord = floor(FirstEscape.Position(1,:));
mask(coord(2)-4:coord(2)+4,coord(1)-4:coord(1)+4)=1;
escapehole = activecontour(I,mask);
Chole = regionprops(escapehole,'Centroid','MajorAxisLength','Area');
delete(FirstEscape);
B = bwboundaries(escapehole);
for i=1:length(B)
    plot(B{i}(:,2),B{i}(:,1),'g','LineWidth',2);
end
slope = (Cmaze.Centroid(2) - Chole.Centroid(2)) ./ (Cmaze.Centroid(1) - Chole.Centroid(1));
if Cmaze.Centroid(1)<Chole.Centroid(1)
angle = atand(slope);
else
angle = 180+atand(slope);    
end
Proximity=0;
Timeintargethole=0;
targetlocationx=Chole.Centroid(1);
targetlocationy=Chole.Centroid(2);
hole = drawcircle('Center',[targetlocationx,targetlocationy],... 
   'Radius',Chole.MajorAxisLength*1.25,'StripeColor','blue');
%%usefull for nsf or something like that
sx=xt(1:(timesec*fps));
sy=yt(1:(timesec*fps));
f=length(sx);
Proximity=[];
for i=1:f-1
d = sqrt((sx(i+1)-targetlocationx)^2+(sy(i+1)-targetlocationy)^2);
Proximity=[Proximity;d];
end
MeanProximity=mean(Proximity*conversionfactor);
Uncertainty=std(Proximity*conversionfactor);
Escaperoi = inROI(hole,sx,sy);
Entries=find(diff(Escaperoi)==1);
No_Entries=length(Entries);
if No_Entries==0
FirstEscape=timesec
else
FirstEscape = min(Entries)/fps
end
DistanceTravelled=sum(Frames)*conversionfactor;
DistanceTravelledtoEscape=sum(Frames(1:(find(Escaperoi, 1, 'first'))))*conversionfactor
Speed=DT;
AverageSpeed=mean(Speed(Speed>3));
if sum(Speed(Speed>9))==0
Rapid_exploratory_behaviour=0;
else
Rapid_exploratory_behaviour=mean(Speed(Speed>9));
end
NosePoking=length(find(diff(Escaperoi)==1));
Timeintargethole=sum(Escaperoi)/fps
Switch=diff(Escaperoi);
% quadrant 
hold on
th = linspace( -45, 45, 10);
Rad = (Cmaze.MajorAxisLength-4*Chole.MajorAxisLength)/2;  %inside according to hole size
Rad = (Cmaze.MajorAxisLength/4);
x = Rad*cosd(th+angle) + Cmaze.Centroid(1);
y = Rad*sind(th+angle) + Cmaze.Centroid(2);
hold on
xCenter = Cmaze.Centroid(1);  
yCenter = Cmaze.Centroid(2);
x = [xCenter, x, xCenter];
y = [yCenter, y, yCenter];
thisROI=[x;y]';
Quadrant=drawpolygon('Position',thisROI,'FaceAlpha',.20,'Color','r');
Qroi = inROI(Quadrant,sx,sy);
TimeintargetQC=sum(Qroi)/fps;
Rad = Cmaze.MajorAxisLength/2;  %or whatever radius you want
x = Rad*cosd(th+angle) + Cmaze.Centroid(1);
y = Rad*sind(th+angle) + Cmaze.Centroid(2);
hold on
xCenter = Cmaze.Centroid(1);  % Let's have the center/tip of the sector be at the middle of the image.
yCenter = Cmaze.Centroid(2);
x = [xCenter, x, xCenter];
y = [yCenter, y, yCenter];
thisROI=[x;y]';
Quadrant=drawpolygon('Position',thisROI,'FaceAlpha',.20,'Color','g');
Qroi = inROI(Quadrant,sx,sy);
TimeintargetQ=(sum(Qroi)/fps)-TimeintargetQC
hold on
plot(filteredcom(1:r,1),filteredcom(1:r,2),'k-','LineWidth',2,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor','none',...
    'MarkerSize',10)
AllData=[DistanceTravelled,DistanceTravelledtoEscape,AverageSpeed,...
    Rapid_exploratory_behaviour,FirstEscape,NosePoking,...
    Timeintargethole,TimeintargetQ,MeanProximity];
%% plotting tracks
imshow(BackgroundCropped)
hold on
plot(x, y, 'g-', 'LineWidth', 2);
hold on
plot(filteredcom(1:r,1),filteredcom(1:r,2),'k-','LineWidth',2,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor','none',...
    'MarkerSize',2)
hold on 
plot(Chole.Centroid(1),Chole.Centroid(2),'g.')
hold on
plot(Cmaze.Centroid(1),Cmaze.Centroid(2),'ro')
hold off
f = gcf;
exportgraphics(f,[filename(1:end-4), '.png'],'Resolution',300)
close all