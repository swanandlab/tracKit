%% Clear previous data
clc; close; clear
FileName=[]
AllData=[]
%% Convert to real world scale(cm)
filename='Trackit MS4_1_2_13_9_19_OFT'
load([filename,'.mat'])
fps=15;
peri=11.375;
AnalyseMin=1;
realworldvalue=45.5;
walking_speed=3;
pixelvalue=sqrt(sum(ClearOutsideMaskCropped,'all'));
conversionfactor=realworldvalue/pixelvalue;
pixperi=floor(peri/conversionfactor);
se = strel('square',2*pixperi);
ClearOutsideMaskCropped = padarray(ClearOutsideMaskCropped,[1 1],0,'both');
inner = imerode(ClearOutsideMaskCropped,se);
imshow(BackgroundCropped)
[Bouter,~,~] = bwboundaries(ClearOutsideMaskCropped);
hold on
boundaryouter = Bouter{1,1};
plot(boundaryouter(:,2), boundaryouter(:,1), 'g','LineWidth',2);
hold on
[Binner,~,~] = bwboundaries(inner);
boundaryinner = Binner{1,1};
plot(boundaryinner(:,2), boundaryinner(:,1), 'r','LineWidth',2);
% Distance Travelled
timesec=AnalyseMin*60;
r =timesec*fps;% for reshaping by frame rate 
dstnc=[];
x=filteredcom(:,1);
y=filteredcom(:,2);
f=length(filteredcom);
for i=1:f-1
d = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2); % calculate distance
dstnc=[dstnc;d];
end
R=dstnc(1:r);
DT = (sum(reshape(R,fps,[]),1))*conversionfactor;
x=Binner{1, 1}(:,1);
y=Binner{1, 1}(:,end);
polyin = polyshape(y,x);
[cx,cy] = centroid(polyin);

plot(filteredcom(1:r,1),filteredcom(1:r,2),'c-','LineWidth',1)
f=gcf;
%saveas(f,[filename, '.tif'])
exportgraphics(f,[filename, '.png'],'Resolution',300)
% Proximity to centre
Proximity=0;
timeintarget=0;
imshow(BackgroundCropped);
targetlocationx=cx;
targetlocationy=cy;
x=filteredcom(:,1);
y=filteredcom(:,2);
sx=x(1:(timesec*fps));
sy=y(1:(timesec*fps));
f=length(sx);
Proximity=[];
for i=1:f-1
d = sqrt((x(i+1)-targetlocationx)^2+(y(i+1)-targetlocationy)^2);
Proximity=[Proximity;d];
end
h=drawpolygon('Position',Binner{1, 1},'FaceAlpha',0,'Color','r');
tf = inROI(h,sx,sy);
close
DistanceTravelled=sum(DT);
DistanceTravelledCentre=sum(R(tf==1))*conversionfactor;
Speed=DT;
AverageSpeed=mean(Speed(Speed>walking_speed));
if sum(Speed(Speed>9))==0
Rapid_exploratory_behaviour=0;
else
Rapid_exploratory_behaviour=mean(Speed(Speed>9));
end
Entries=find(diff(tf)==1);
No_Entries=length(Entries);
if No_Entries==0
FirstEntry=timesec;
else
FirstEntry = min(Entries)/fps;
end
%if mouse put in centre then first reentry
Exits=find(diff(tf)==-1);
No_Exits=length(Exits);
if No_Exits==0
FirstExit=0;
else
FirstExit = min(Exits)/fps;
end
timeintarget=sum(tf)/fps;
Switch=diff(tf);
MeanProximity=mean(Proximity*conversionfactor);
FileName=[FileName;{filename}]
AllData=[AllData;([DistanceTravelled,DistanceTravelledCentre,AverageSpeed,...
    Rapid_exploratory_behaviour,FirstEntry,FirstExit,No_Entries,...
    No_Exits,timeintarget,MeanProximity])]