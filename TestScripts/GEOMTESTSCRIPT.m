clear all
close all
clc

%% Test with shapes that intersect at points

%% Test with shapes that interscect on curves
g1 = Geometry2D.draw('Sector','Radius',[1 2],'Angle',2*pi*7/8);
g2 = Geometry2D.draw('Sector','Radius',[1 2],'Angle',2*pi/8,'Rotation',2*pi*7/8);
g3 = Geometry2D.draw('Sector','Radius',[1 2],'Angle',2*pi/8,'Rotation',-2*pi/8);

assert(abs(g1.area - pi*7/8*3) < sqrt(eps));
assert(abs(g2.area - pi/8*3)   < sqrt(eps));
assert(abs(g3.area - pi/8*3)   < sqrt(eps));

s  = g1+g2;
assert(abs(s.area - pi*3) < sqrt(eps));

s  = g1+g3;
assert(abs(s.area - pi*3) < sqrt(eps));

g4 = Geometry2D.draw('Rectangle','Length',0.5,'Width',0.5,'Base','Center','Position',[1.5,-0.25]);

s  = g1+g4;
assert(abs(s.area - pi*7/8*3-1/4) < sqrt(eps));

g5 = Geometry2D.draw('Rectangle','Length',0.5,'Width',0.5,'Base','Center','Position',[1.5,+0.25]);

s  = g2+g5;
assert(abs(s.area - pi/8*3-1/4) < sqrt(eps));

s  = g3+g5;
assert(abs(s.area - pi/8*3-1/4) < sqrt(eps));

s  = g4+g5;
assert(abs(s.area - 1/2) < sqrt(eps));

g6 = Geometry2D.draw('Sector','Radius',[2 3],'Angle',2*pi*7/8);

s  = g1+g6;
assert(abs(s.area - pi*7/8*8) < sqrt(eps));

g7 = Geometry2D.draw('Sector','Radius',[2 3],'Angle',2*pi*5/8,'Rotation',2*pi/8);

s  = g1+g7;
assert(abs(s.area - pi*7/8*3 - pi*5/8*5)<sqrt(eps));

g8 = g5.rotate('Position',[0 0],'Rotation',-2*pi/8);

s  = g1+g8;
assert(abs(s.area - pi*7/8*3-1/4) < sqrt(eps));

%% Test with shapes in which one fully overlaps the other
g1 = [Geometry2D.draw('Rectangle','Length',1,'Width',1),...
      Geometry2D.draw('Sector','Radius',[1 3],'Angle',2*pi,'Position',[0 -2])];

g2 = [Geometry2D.draw('Rectangle','Length',0.5,'Width',0.5),...
      Geometry2D.draw('Sector','Radius',[0.125,0.25],'Angle',pi/8)];
  
assert(abs(g1(1).area - 1)    < sqrt(eps));
assert(abs(g1(2).area - pi*8) < sqrt(eps));
assert(abs(g2(1).area - 0.25) < sqrt(eps));
assert(abs(g2(2).area - pi/16*(0.25^2-0.125^2)) < sqrt(eps));

d        = zeros(2,2,3);
d(:,:,1) = [1 2;
            1 2];
d(:,:,2) = [1 1;
            1 1];
d(:,:,3) = [2 3;
            2 3];

c        = zeros(2,2,3);
c(:,:,1) = [4 2;
            4 2];
c(:,:,2) = [4 4;
            4 4];
c(:,:,3) = [8 6;
            8 6];

for j = 1:numel(g1)
    for i = 1:numel(g2);
        [j i];
        s1 = g1(j)+g2(i);
        assert(numel(s1.Domains) == d(i,j,1));
        assert(length(s1.Curves) == c(i,j,1));
        assert(abs(s1.area-g1(j).area) < sqrt(eps));

        i1 = g1(j)*g2(i);
        assert(numel(i1.Domains) == d(i,j,2));
        assert(length(i1.Curves) == c(i,j,2));
        assert(abs(i1.area-g2(i).area) < sqrt(eps));

        d1 = g1(j)-g2(i);
        assert(numel(d1.Domains) == d(i,j,3));
        assert(length(d1.Curves) == c(i,j,3));
        assert(abs(d1.area-g1(j).area+g2(i).area) < sqrt(eps));

        i2 = d1*g2(i);
        assert(isempty(i2.Domains{1}));
        assert(isempty(i2.Curves));
        assert(abs(i2.area) < sqrt(eps));

        d2 = g2(i)-g1(j);
        assert(isempty(d2.Domains{1}));
        assert(isempty(d2.Curves));
        assert(abs(d2.area) < sqrt(eps));
    end
end

%% Test with partially overlaping shapes

% clear all
% close all
% clc
% 
% tic
% S = MotorProto('TEST');
% toc
% 
% tic
% P = S.Parameters;
% P.new('np',8);
% P.new('nph',3);
% P.new('nspp',1);
% P.new('swp',0.5);
% P.new('tf1',0.9);
% P.new('bil',0.4);
% P.new('sll',0.55);
% P.new('rso',12);
% P.new('slfrac',0.6);
% P.new('lag',0.001);
% P.new('thetat1',2*pi/np/nph/nspp);
% P.new('thetat2',2*pi/np/nph/nspp*tf1/2);
% P.new('thetat3',2*pi/np/nph/nspp*(1-tf1/2));
% P.new('rsi',rso*slfrac);
% P.new('swm',2*rsi*sin(thetat1*swp));
% P.new('sm',-1/tan(thetat1/2));
% P.new('sdx1',swm/sqrt(1+sm^2)*(1-tf1)/2);
% P.new('sdx2',swm/sqrt(1+sm^2)*swp/2);
% P.new('sdy1',sm*sdx1);
% P.new('sdy2',sm*sdx2);
% P.new('rtf',rsi*(bil+sll)+rso*(1-bil-sll));
% P.new('rbi',rsi*bil+rso*(1-bil));
% P.new('rro',rsi-lag);
% P.new('rri',rro*0.5);
% P.new('xo',0);
% P.new('yo',0);
% toc
% 
% tic
% P.edit('np',10);
% P.edit('np',nph*2);
% toc
% 
% tic
% Geometry0D.draw('point','X',rso,'Y',0);
% Geometry0D.draw('point','X',rso*cos(thetat1),'Y',rso*sin(thetat1));
% Geometry0D.draw('point','X',rsi,'Y',0);
% Geometry0D.draw('point','X',rsi*cos(thetat1),'Y',rsi*sin(thetat1));
% Geometry0D.draw('point','X',rro,'Y',0);
% Geometry0D.draw('point','X',rro*cos(2*pi/np),'Y',rro*sin(2*pi/np));
% Geometry0D.draw('point','X',rri,'Y',0);
% Geometry0D.draw('point','X',rri*cos(2*pi/np),'Y',rri*sin(2*pi/np));
% toc
% 
% tic
% P.edit('np',8);
% P.edit('nspp',2);
% toc
%  
% tic
% statorSlot = [];
% statorToothGap = [];
% statorIron= [];
% rotorIron = [];
% 
% % %test line drawing
% C{1} = Geometry1D.draw('Arc','Name','c1','Position',[0 0],'Radius',rso,'Angle',thetat1,'Rotation',0);
% 	statorIron = [statorIron,1];
% C{2} = Geometry1D.draw('Arc','Position',[0 0],'Radius',rsi,'Angle',thetat2,'Rotation',thetat3);
%      statorIron = [statorIron,2];
% C{3} = Geometry1D.draw('Arc','Position',[0 0],'Radius',rsi,'Angle',thetat2,'Rotation',0);
%      statorIron = [statorIron,3];
% C{4} = Geometry1D.draw('Arc','Position',[0 0],'Radius',rsi,'Angle',thetat3-thetat2,'Rotation',thetat2);
%      statorToothGap = [statorToothGap,4];
% C{5} = Geometry1D.draw('Polyline','Points',[rso*cos(thetat1) rso*sin(thetat1);rsi*cos(thetat1) rsi*sin(thetat1)]);
%      statorIron = [statorIron,5];
% C{6} = Geometry1D.draw('Polyline','Points',[rsi 0;rso 0]);
%     statorIron = [statorIron,6];
% C{7} = Geometry1D.draw('Polyline','Points',[rtf*cos(thetat1/2)+sdx1 rtf*sin(thetat1/2)+sdy1;rtf*cos(thetat1/2)-sdx1 rtf*sin(thetat1/2)-sdy1]);
%     statorToothGap = [statorToothGap,7];
%     statorSlot = [statorSlot,7];
% C{8} = Geometry1D.draw('Polyline','Points',[rtf*cos(thetat1/2)+sdx1 rtf*sin(thetat1/2)+sdy1;rtf*cos(thetat1/2)+sdx2 rtf*sin(thetat1/2)+sdy2]);
%     statorIron = [statorIron,8];
%     statorSlot = [statorSlot,8];
% C{9} = Geometry1D.draw('Polyline','Points',[rtf*cos(thetat1/2)-sdx1 rtf*sin(thetat1/2)-sdy1;rtf*cos(thetat1/2)-sdx2 rtf*sin(thetat1/2)-sdy2]);
%     statorIron = [statorIron,9];
%     statorSlot = [statorSlot,9];
% C{10} = Geometry1D.draw('Polyline','Points',[rtf*cos(thetat1/2)-sdx1 rtf*sin(thetat1/2)-sdy1;rsi*cos(thetat3) rsi*sin(thetat3)]);
%     statorToothGap = [statorToothGap,10];
%     statorIron = [statorIron,10];
% C{11} = Geometry1D.draw('Polyline','Points',[rtf*cos(thetat1/2)+sdx1 rtf*sin(thetat1/2)+sdy1;rsi*cos(thetat2) rsi*sin(thetat2)]);
%     statorIron = [statorIron,11];
%     statorToothGap = [statorToothGap,11];
% C{12} = Geometry1D.draw('Polyline','Points',[rbi*cos(thetat1/2)+sdx2 rbi*sin(thetat1/2)+sdy2;rbi*cos(thetat1/2)-sdx2 rbi*sin(thetat1/2)-sdy2]);
%     statorIron = [statorIron,12];
%     statorSlot = [statorSlot,12];
% C{13} = Geometry1D.draw('Polyline','Points',[rbi*cos(thetat1/2)+sdx2 rbi*sin(thetat1/2)+sdy2;rtf*cos(thetat1/2)+sdx2 rtf*sin(thetat1/2)+sdy2]);
%     statorIron = [statorIron,13];
%     statorSlot = [statorSlot,13];
% C{14} = Geometry1D.draw('Polyline','Points',[rbi*cos(thetat1/2)-sdx2 rbi*sin(thetat1/2)-sdy2;rtf*cos(thetat1/2)-sdx2 rtf*sin(thetat1/2)-sdy2]);
%     statorIron = [statorIron,14];
%     statorSlot = [statorSlot,14];
% C{15} = Geometry1D.draw('Arc','Position',[0 0]','Radius',rro,'Angle',2*pi/np,'Rotation',0);
%     rotorIron = [rotorIron,15];
% C{16} = Geometry1D.draw('Polyline','Points',[rro*cos(2*pi/np) rro*sin(2*pi/np);rri*cos(2*pi/np) rri*sin(2*pi/np)]);
%     rotorIron = [rotorIron,16];
% C{17} = Geometry1D.draw('Arc','Position',[0 0],'Radius',rri,'Angle',2*pi/np,'Rotation',0);
%     rotorIron = [rotorIron,17];
% C{18} = Geometry1D.draw('Polyline','Points',[rri 0;rro 0]);
%     rotorIron = [rotorIron,18];
% toc
% 
% tic
% P.edit('np',6);
% P.edit('nspp',3);
% toc
% 
% %test some intersections methods
% tic
% C{19}=Geometry1D.draw('Polyline','Points',[rri*0.0 0;rri*1.5*cos(pi/np) rri*1.5*sin(pi/np)]);
% C{20}=Geometry1D.draw('Arc','Position',[rri -rri/2],'Radius',rri,'Angle',2*pi/np,'Rotation',2*pi/np);
% C{21}=Geometry1D.draw('Polyline','Points',[rri*1.05 -rri*0.05;rri*1.05*cos(2*pi/np) rri*1.05*sin(2*pi/np)]);
% C{22}=Geometry1D.draw('Arc','Position',[rri+rri/2*1.1 rri],'Radius','rri','Angle','2*pi/np','Rotation',pi+pi/np/4);
% C{23}=Geometry1D.draw('Arc','Position',[0 0],'Radius',rri,'Angle',2*pi/np,'Rotation',-0.1);
% C{24}=Geometry1D.draw('Arc','Position',[0 0],'Radius',rri,'Angle',pi/np,'Rotation',0.1);
% C{25}=Geometry1D.draw('Arc','Position',[0 0],'Radius',rri,'Angle',2*pi/np,'Rotation',2*pi/np);
% C{26}=Geometry1D.draw('Polyline','Points',[0 0;rri rro]);
% C{27}=Geometry1D.draw('Polyline','Points',[rri rro;0 0]);
% C{28}=Geometry1D.draw('Polyline','Points',[0 0;rri/2 rro/2]);
% C{29}=Geometry1D.draw('Polyline','Points',[rri/4 rro/4;rri*3/4 rro*3/4]);
% C{30}=Geometry1D.draw('Polyline','Points',[0 0;0 rro]);
% toc
% 
% tic
% S=intersection(C{17},C{18});%arc line endpoint
% S=intersection(C{17},C{19});%arc line one point
% S=intersection(C{17},C{20});%arc arc one point
% S=intersection(C{17},C{21});%arc line two points
% S=intersection(C{17},C{22});%arc arc two points
% S=intersection(C{19},C{21});%line line intersection
% S=intersection(C{17},C{17});%indentical arcs
% S=intersection(C{17},C{23});%partially overlapping arcs
% S=intersection(C{17},C{24});%subset arcs
% S=intersection(C{17},C{25});%arc arc endpoint
% S=intersection(C{26},C{26});%identical lines
% S=intersection(C{26},C{27});%identical lines, opposite directions
% S=intersection(C{26},C{28});%subset lines with endpoints
% S=intersection(C{26},C{29});%subset without endpoints
% S=intersection(C{28},C{29});%overlapping lines
% S=intersection(C{28},C{30});%line line endpoint
% toc
% 
% tic
% P.edit('np',6);
% P.edit('nspp',3);
% toc
% 
% %create planes from curves
% % tic
% % Geometry2D.draw('Curvewise','Curves',{C{statorIron}},'PlotStyle',{'b'});
% % Geometry2D.draw('Curvewise','Curves',{C{statorSlot}},'PlotStyle',{'y'});
% % Geometry2D.draw('Curvewise','Curves',{C{statorToothGap}},'PlotStyle',{'w'});
% % Geometry2D.draw('Curvewise','Curves',{C{rotorIron}},'PlotStyle',{'b'});
% % toc
% 
% tic
% %draw some rectangles
% r1=Geometry2D.draw('rect','Length',rri,'Width',rsi,'Rotation',0,'Position',[-rri 0]);
% r2=Geometry2D.draw('rect','Length',rri,'Width',rsi,'Rotation',0,'Position',[-rri/4 rsi/4],'Base','Center');
% toc
% 
% tic
% %draw some sectors
% s1=Geometry2D.draw('sector','Radius',[rri rro],'Position',[-rsi 0],'Angle',pi);
% s2=Geometry2D.draw('sector','Radius',[rri rro],'Position',[-rsi*2/3 0],'Angle',pi,'Rotation',-pi/2);
% toc
% 
% tic
% P.edit('rri',1);
% P.edit('rri',4);
% toc
% 
% %test unions
% tic
% u1 = r1+r2;clf;u1.plot;pause
% assert(numel(u1.Domains) == 1);
% assert(numel(u1.Curves)  == 8);
% 
% u2 = r1+s1;clf;u2.plot;pause
% assert(numel(u1.Domains) == 1);
% assert(numel(u1.Curves)  == 7);
% 
% u3 = r1+s2;clf;u3.plot;pause
% assert(numel(u1.Domains) == 2);
% assert(numel(u1.Curves)  == 8);
% 
% u4 = r2+s1;clf;u4.plot;pause
% u5 = r2+s2;clf;u5.plot;pause
% u6 = s1+s2;clf;u6.plot;pause
% toc
% 
% %test intersections
% 
% tic
% i1 = r1*r2;clf;i1.plot;pause
% i2 = r1*s1;clf;i2.plot;pause
% i3 = r1*s2;clf;i3.plot;pause
% i4 = r2*s1;clf;i4.plot;pause
% i5 = r2*s2;clf;i5.plot;pause
% i6 = s1*s2;clf;i6.plot;pause
% toc
% 
% 
% %test differences
% tic
% d1 = r1-r2;clf;d1.plot;pause
% d2 = r1-s1;clf;d2.plot;pause
% d3 = r1-s2;clf;d3.plot;pause
% d4 = r2-r1;clf;d4.plot;pause
% d5 = r2-s1;clf;d5.plot;pause
% d6 = r2-s2;clf;d6.plot;pause
% d7 = s1-r1;clf;d7.plot;pause
% d8 = s1-r2;clf;d8.plot;pause
% d9 = s1-s2;clf;d9.plot;pause
% d10 = s2-s1;clf;d10.plot;pause
% toc
% 
% a = {r1,r2,s1,s2};
% D = r1+r2+s1+s2;
% for i=1:4
%     t=D+a{i};clf;t.plot;pause
%     t=D-a{i};clf;t.plot;pause
%     t=D*a{i};clf;t.plot;pause
%     t=(D-a{i})+a{i};clf;t.plot;pause
% end
% 
% % %now draw the same geometry using plane methods
% tic
% iron=G.draw('sector','Radius','[rsi,rso]','Angle','thetat1','Position','[0 0]','PlotStyle',{'b'},'Rotation','thetat1');
% slot=G.draw('rect','Length','(rso-rsi)*sll','Width','swm/2','Rotation','thetat1*3/2','Position','[(rtf+(rso-rsi)*sll/2)*cos(thetat1*3/2),(rtf+(rso-rsi)*sll/2)*sin(thetat1*3/2)]','base','center','PlotStyle',{'y'});
% gap=G.draw('rect','Length','(rso-rsi)*sll','Width','swm/2*.2','Rotation','thetat1*3/2','Position','[(rtf+(rso-rsi)*sll/2)*cos(thetat1*3/2)*.9,(rtf+(rso-rsi)*sll/2)*sin(thetat1*3/2)*.9]','base','center','PlotStyle',{'w'});
% iron = iron-slot;
% gap = gap*iron;
% iron = iron-gap;
% iron.plot;
% gap.plot;
% slot.plot;