% Plot Synthetic Seismogram with Matlab
% 
% Daniel Koehn
% Raisdorf, den 11.3.2005

clear all
close all

%load('full_wave_test_start_40_x.mat')
%sek1=sek;

FW=0.0;

%t1=1:ns;
%t1=t1.*dt;


%amp=1e-1;

% Calculate time scale
%t=1:ns;
%t=t.*dt;

%XREC1=170.0;
%YREC1=80.0;

%XREC2=170.0;
%YREC2=270.0;

XREC1=215.0;
YREC1=155.0;

XREC2=4785.0;
YREC2=155.0;

XSOUR1=215.0;
YSOUR1=20.0;

XSOUR2=4785.0;
YSOUR2=20.0;

%XREC1=130.0+FW;
%YREC1=0.0;

%XREC2=130.0+FW;
%YREC2=360.0;

drec=15.0;
dsour=225.0;

% receiver
ntr=301;

for i=1:ntr

% Zeitrichtung umkehren
h=0;
%for j=1:ns
%sek2(j,i)=sek1(ns-h,i);
%h=h+1;
%end

% Berechne Geophonkoordinaten
YREC(i)=YREC1;
XREC(i)=XREC1+(i-1).*drec;
ZREC(i)=0.0;

TD(i)=0.0;
FMAIN(i)=30.0;
AMP(i)=1.0;

end

%geoph=[XREC' ZREC' YREC' TD' FMAIN' AMP'];

%dlmwrite('source.dat',geoph,'delimiter','\t','precision','%.2f');

geoph1=[XREC' YREC'];

dlmwrite('receiver_resize.dat',geoph1,'delimiter','\t','precision','%.2f');


% sources
ntr=21;

for i=1:ntr

% Berechne Geophonkoordinaten
YSOUR(i)=YSOUR1;
XSOUR(i)=XSOUR1+(i-1).*dsour;
ZSOUR(i)=0.0;

TD(i)=0.0;
FMAIN(i)=20.0;
AMP(i)=1.0;


pause(0.1);

end

% Ausgabe der Quellpunktkoordinaten

%for i=1:ntr
%imfile=['source_resize_',int2str(i),'.dat'];
%fid = fopen(imfile,'w');
%fprintf(fid,'%i \n',t1);
%fprintf(fid,'%6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \n',XSOUR(i),ZSOUR(i),YSOUR(i),TD(i),FMAIN(i),AMP(i));
%fclose(fid);
%end

imfile=['sources_resize.dat'];
fid = fopen(imfile,'w');
for i=1:ntr
fprintf(fid,'%6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f \n',XSOUR(i),ZSOUR(i),YSOUR(i),TD(i),FMAIN(i),AMP(i));
end

fclose(fid);
