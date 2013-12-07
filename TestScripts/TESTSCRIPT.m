clear all
close all
clc

tic
S = SSSolve('TEST');
toc

tic
P = S.parameters;
P.new('np',8);
P.new('nph',3);
P.new('nspp',1);
P.new('swp',0.5);
P.new('tf1',0.9);
P.new('bil',0.4);
P.new('sll',0.55);
P.new('rso',12);
P.new('slfrac',0.6);
P.new('lag',0.001);
P.new('thetat1',2*pi/np/nph/nspp);
P.new('thetat2',2*pi/np/nph/nspp*tf1/2);
P.new('thetat3',2*pi/np/nph/nspp*(1-tf1/2));
P.new('rsi',rso*slfrac);
P.new('swm',2*rsi*sin(thetat1*swp));
P.new('sm',-1/tan(thetat1/2));
P.new('sdx1',swm/sqrt(1+sm^2)*(1-tf1)/2);
P.new('sdx2',swm/sqrt(1+sm^2)*swp/2);
P.new('sdy1',sm*sdx1);
P.new('sdy2',sm*sdx2);
P.new('rtf',rsi*(bil+sll)+rso*(1-bil-sll));
P.new('rbi',rsi*bil+rso*(1-bil));
P.new('rro',rsi-lag);
P.new('rri',rro*0.5);
P.new('xo',0);
P.new('yo',0);
toc

tic
P.edit('np',10);
P.edit('np',nph*2);
toc

tic
M = S.models;
G = M.geometry;
G.rtPlot = true;
G.draw('point','x',rso,'y',0);
G.draw('point','x',rso*cos(thetat1),'y',rso*sin(thetat1));
G.draw('point','x',rsi,'y',0);
G.draw('point','x',rsi*cos(thetat1),'y',rsi*sin(thetat1));
G.draw('point','x',rro,'y',0);
G.draw('point','x',rro*cos(2*pi/np),'y',rro*sin(2*pi/np));
G.draw('point','x',rri,'y',0);
G.draw('point','x',rri*cos(2*pi/np),'y',rri*sin(2*pi/np));
toc

tic
P.edit('np',8);
P.edit('nspp',2);
toc
 
tic
statorSlot = [];
statorToothGap = [];
statorIron= [];
rotorIron = [];

% %test line drawing
G.draw('arc','name','c1','position',[0 0],'radius',rso,'angle',thetat1,'rotation',0);
	statorIron = [statorIron,1];
G.draw('arc','position',[0 0],'radius',rsi,'angle',thetat2,'rotation',thetat3);
     statorIron = [statorIron,2];
G.draw('arc','position',[0 0],'radius',rsi,'angle',thetat2,'rotation',0);
     statorIron = [statorIron,3];
G.draw('arc','position',[0 0],'radius',rsi,'angle',thetat3-thetat2,'rotation',thetat2);
     statorToothGap = [statorToothGap,4];
G.draw('polyline','points',[rso*cos(thetat1) rso*sin(thetat1);rsi*cos(thetat1) rsi*sin(thetat1)]);
     statorIron = [statorIron,5];
G.draw('polyline','points',[rsi 0;rso 0]);
    statorIron = [statorIron,6];
G.draw('polyline','points',[rtf*cos(thetat1/2)+sdx1 rtf*sin(thetat1/2)+sdy1;rtf*cos(thetat1/2)-sdx1 rtf*sin(thetat1/2)-sdy1]);
    statorToothGap = [statorToothGap,7];
    statorSlot = [statorSlot,7];
G.draw('polyline','points',[rtf*cos(thetat1/2)+sdx1 rtf*sin(thetat1/2)+sdy1;rtf*cos(thetat1/2)+sdx2 rtf*sin(thetat1/2)+sdy2]);
    statorIron = [statorIron,8];
    statorSlot = [statorSlot,8];
G.draw('polyline','points',[rtf*cos(thetat1/2)-sdx1 rtf*sin(thetat1/2)-sdy1;rtf*cos(thetat1/2)-sdx2 rtf*sin(thetat1/2)-sdy2]);
    statorIron = [statorIron,9];
    statorSlot = [statorSlot,9];
G.draw('polyline','points',[rtf*cos(thetat1/2)-sdx1 rtf*sin(thetat1/2)-sdy1;rsi*cos(thetat3) rsi*sin(thetat3)]);
    statorToothGap = [statorToothGap,10];
    statorIron = [statorIron,10];
G.draw('polyline','points',[rtf*cos(thetat1/2)+sdx1 rtf*sin(thetat1/2)+sdy1;rsi*cos(thetat2) rsi*sin(thetat2)]);
    statorIron = [statorIron,11];
    statorToothGap = [statorToothGap,11];
G.draw('polyline','points',[rbi*cos(thetat1/2)+sdx2 rbi*sin(thetat1/2)+sdy2;rbi*cos(thetat1/2)-sdx2 rbi*sin(thetat1/2)-sdy2]);
    statorIron = [statorIron,12];
    statorSlot = [statorSlot,12];
G.draw('polyline','points',[rbi*cos(thetat1/2)+sdx2 rbi*sin(thetat1/2)+sdy2;rtf*cos(thetat1/2)+sdx2 rtf*sin(thetat1/2)+sdy2]);
    statorIron = [statorIron,13];
    statorSlot = [statorSlot,13];
G.draw('polyline','points',[rbi*cos(thetat1/2)-sdx2 rbi*sin(thetat1/2)-sdy2;rtf*cos(thetat1/2)-sdx2 rtf*sin(thetat1/2)-sdy2]);
    statorIron = [statorIron,14];
    statorSlot = [statorSlot,14];
G.draw('arc','position',[0 0]','radius',rro,'angle',2*pi/np,'rotation',0);
    rotorIron = [rotorIron,15];
G.draw('polyline','points',[rro*cos(2*pi/np) rro*sin(2*pi/np);rri*cos(2*pi/np) rri*sin(2*pi/np)]);
    rotorIron = [rotorIron,16];
G17=G.draw('arc','position',[0 0],'radius',rri,'angle',2*pi/np,'rotation',0);
    rotorIron = [rotorIron,17];
G18=G.draw('polyline','points',[rri 0;rro 0]);
    rotorIron = [rotorIron,18];
toc

tic
P.edit('np',6);
P.edit('nspp',3);
toc

%test some intersections methods
tic
G19=G.draw('polyline','points',[rri*0.0 0;rri*1.5*cos(pi/np) rri*1.5*sin(pi/np)]);
G20=G.draw('arc','position',[rri -rri/2],'radius',rri,'angle',2*pi/np,'rotation',2*pi/np);
G21=G.draw('polyline','points',[rri*1.05 -rri*0.05;rri*1.05*cos(2*pi/np) rri*1.05*sin(2*pi/np)]);
G22=G.draw('arc','position',[rri+rri/2*1.1 rri],'radius','rri','angle','2*pi/np','rotation',pi+pi/np/4);
G23=G.draw('arc','position',[0 0],'radius',rri,'angle',2*pi/np,'rotation',-0.1);
G24=G.draw('arc','position',[0 0],'radius',rri,'angle',pi/np,'rotation',0.1);
G25=G.draw('arc','position',[0 0],'radius',rri,'angle',2*pi/np,'rotation',2*pi/np);
G26=G.draw('polyline','points',[0 0;rri rro]);
G27=G.draw('polyline','points',[rri rro;0 0]);
G28=G.draw('polyline','points',[0 0;rri/2 rro/2]);
G29=G.draw('polyline','points',[rri/4 rro/4;rri*3/4 rro*3/4]);
G30=G.draw('polyline','points',[0 0;0 rro]);
toc

tic
S=intersection(G.curves{17},G.curves{18});%arc line endpoint
S=intersection(G.curves{17},G.curves{19});%arc line one point
S=intersection(G.curves{17},G.curves{20});%arc arc one point
S=intersection(G.curves{17},G.curves{21});%arc line two points
S=intersection(G.curves{17},G.curves{22});%arc arc two points
S=intersection(G.curves{19},G.curves{21});%line line intersection
S=intersection(G.curves{17},G.curves{17});%indentical arcs
S=intersection(G.curves{17},G.curves{23});%partially overlapping arcs
S=intersection(G.curves{17},G.curves{24});%subset arcs
S=intersection(G.curves{17},G.curves{25});%arc arc endpoint
S=intersection(G.curves{26},G.curves{26});%identical lines
S=intersection(G.curves{26},G.curves{27});%identical lines, opposite directions
S=intersection(G.curves{26},G.curves{28});%subset lines with endpoints
S=intersection(G.curves{26},G.curves{29});%subset without endpoints
S=intersection(G.curves{28},G.curves{29});%overlapping lines
S=intersection(G.curves{28},G.curves{30});%line line endpoint
toc

tic
P.edit('np',6);
P.edit('nspp',3);
toc

%create planes from curves
tic
G.draw('curvewise','curves',{G.curves{statorIron}},'plotStyle',{'b'});
G.draw('curvewise','curves',{G.curves{statorSlot}},'plotStyle',{'y'});
G.draw('curvewise','curves',{G.curves{statorToothGap}},'plotStyle',{'w'});
G.draw('curvewise','curves',{G.curves{rotorIron}},'plotStyle',{'b'});
toc

tic
%draw some rectangles
r1=G.draw('rect','length',rri,'width',rsi,'rotation',0,'position',[-rri 0]);
r2=G.draw('rect','length',rri,'width',rsi,'rotation',0,'position',[-rri/4 rsi/4],'base','center');
toc

tic
%draw some sectors
s1=G.draw('sector','radius',[rri rro],'position',[-rsi 0],'angle',pi);
s2=G.draw('sector','radius',[rri rro],'position',[-rsi*2/3 0],'angle',pi,'rotation',-pi/2);
toc

tic
P.edit('rri',1);
P.edit('rri',4);
toc

%test unions
tic
u1 = r1+r2;clf;u1.plot;pause
u2 = r1+s1;clf;u2.plot;pause
u3 = r1+s2;clf;u3.plot;pause
u4 = r2+s1;clf;u4.plot;pause
u5 = r2+s2;clf;u5.plot;pause
u6 = s1+s2;clf;u6.plot;pause
toc

%test intersections

tic
i1 = r1*r2;clf;i1.plot;pause
i2 = r1*s1;clf;i2.plot;pause
i3 = r1*s2;clf;i3.plot;pause
i4 = r2*s1;clf;i4.plot;pause
i5 = r2*s2;clf;i5.plot;pause
i6 = s1*s2;clf;i6.plot;pause
toc


%test differences
tic
d1 = r1-r2;clf;d1.plot;pause
d2 = r1-s1;clf;d2.plot;pause
d3 = r1-s2;clf;d3.plot;pause
d4 = r2-r1;clf;d4.plot;pause
d5 = r2-s1;clf;d5.plot;pause
d6 = r2-s2;clf;d6.plot;pause
d7 = s1-r1;clf;d7.plot;pause
d8 = s1-r2;clf;d8.plot;pause
d9 = s1-s2;clf;d9.plot;pause
d10 = s2-s1;clf;d10.plot;pause
toc

a = {r1,r2,s1,s2};
D = r1+r2+s1+s2;
for i=1:4
    t=D+a{i};clf;t.plot;pause
    t=D-a{i};clf;t.plot;pause
    t=D*a{i};clf;t.plot;pause
    t=(D-a{i})+a{i};clf;t.plot;pause
end

% %now draw the same geometry using plane methods
tic
iron=G.draw('sector','radius','[rsi,rso]','angle','thetat1','position','[0 0]','plotStyle',{'b'},'rotation','thetat1');
slot=G.draw('rect','length','(rso-rsi)*sll','width','swm/2','rotation','thetat1*3/2','position','[(rtf+(rso-rsi)*sll/2)*cos(thetat1*3/2),(rtf+(rso-rsi)*sll/2)*sin(thetat1*3/2)]','base','center','plotStyle',{'y'});
gap=G.draw('rect','length','(rso-rsi)*sll','width','swm/2*.2','rotation','thetat1*3/2','position','[(rtf+(rso-rsi)*sll/2)*cos(thetat1*3/2)*.9,(rtf+(rso-rsi)*sll/2)*sin(thetat1*3/2)*.9]','base','center','plotStyle',{'w'});
iron = iron-slot;
gap = gap*iron;
iron = iron-gap;
iron.plot;
gap.plot;
slot.plot;