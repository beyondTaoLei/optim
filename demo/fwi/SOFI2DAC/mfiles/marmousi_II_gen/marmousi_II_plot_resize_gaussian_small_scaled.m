close all
clear all

% This m-file is for making movies of wave propagation.

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%
% Here the basic name of the binary snapshot files must be given:
% (The default extension is *.bin. The increasing number of the snapshot
% is added to the basic name, e.g. if ten snapshots were computed
% and the basic name is snap, then the filenames for the snapshots
% are snap1.bin, snap2.bin, ..., snap10.bin.
%filediv='snap/Khang_hom.bin.p.00';

% Receiver Positions
%yrec=3.6:0.66:69.6;
%xrec=45.5.*ones(length(yrec),1);
clear all
close all

% plot distance of sources and receivers
dnrec=1;
dnsour=1;

zeropad=23;     % size of the padding layer for the Gaussian filter 
water_depth=29; % water depth in grid points
filter_depth=40; % depth in which the Gaussian filter is applied 

% write model output files
WRITEMODE=1;

rec=load('receiver_resize.dat');
xrec=rec(:,1);
yrec=rec(:,2);

xrec = xrec(1:dnrec:length(xrec));
yrec = yrec(1:dnrec:length(yrec));

source=load('sources_plot_resize.dat');
xshot=source(:,1);
yshot=source(:,3);

xshot = xshot(1:dnsour:length(xshot));
yshot = yshot(1:dnsour:length(yshot));

% -------------------------------------------------------------------------
% P-Wave Velocity
% -------------------------------------------------------------------------

figure;

% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=13601;
NY1=1; NY2=2801; 
IDX=1; IDY=1;
dh=1.25;

SMOOTH=1;
SHOW=3;
IMP=1;

% time increment for snapshots:
TSNAPINC=2.0e-2; TSNAP1=0.001;
FW=0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

% clip values for Pressure P 
%caxis_value_p_1=2178.0;
%caxis_value_p=1100.0;

% clip values for vx
%caxis_value_p_1=3000.0;
%caxis_value_p=1000.0;

%caxis_value_p_1=2350.0;
%caxis_value_p=900.0;

       
% source position
%xs=2610.0;
%ys=2100.0;

% receiver position
%xrec(1)=2610.0;
%xrec(2)=3730.0;
%yrec(1)=3230.0;
%yrec(2)=2100.0;

%firstframe=50;
%lastframe=150;

%load 'seismic.map'
colormap(jet(256));

% load model
 file='Model_vp.bin';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 vp=fread(fid,[ny,nx],'float');
 fclose(fid);
 
 % cut water column
 vp_wo_water = vp(246:2800,:);
 
 clear vp;
 
 size(vp_wo_water)
 
 % resample model
 vp_resamp = vp_wo_water(1:2:2475,1:2:13600);
 
 clear vp_wo_water;
 clear x;
 clear y;

vp_resamp_1 = vp_resamp(1:2:1160,2801:2:4800);
clear vp_resamp;
vp_resamp = vp_resamp_1;

size(vp_resamp)

caxis_value_vp1 = min(min(vp_resamp))
caxis_value_vp2 = max(max(vp_resamp))
 
% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=1000;
NY1=1; NY2=580; 
IDX=1; IDY=1;
dh=5.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1 
ny=NY2-NY1+1

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

if(SMOOTH==1)
% apply median filter
% add boundary to the model to avoid filter artefacts on the inversion grid 
ZIs = ones(ny+2.*zeropad,nx+2.*zeropad);

for j=1:ny
ZIs(j+zeropad,1:zeropad) = vp_resamp(j,1);
ZIs(j+zeropad,nx+zeropad:nx+2.*zeropad) = vp_resamp(j,1);
end

ZIs(1+zeropad:ny+zeropad,1+zeropad:nx+zeropad) = vp_resamp(1:ny,1:nx); 

for i=1:nx+(2.*zeropad)
ZIs(1:zeropad,i) = ZIs(1+zeropad,i);
ZIs(ny+zeropad:ny+2.*zeropad,i) = ZIs(ny+zeropad,i);
end

sigma = zeropad;

for j=1:ny+2.*zeropad
    for i=1:nx+2.*zeropad
       ZIss(j,i) = ZIs(j,i);
    end
end

% apply gaussian filter
for j=1+zeropad+filter_depth:ny+zeropad
    for i=1+zeropad:nx+zeropad
        
        ZIss(j,i)=0.0;
        normgauss=0.0;
        
        for j1=j-zeropad:j+zeropad
           for i1=i-zeropad:i+zeropad
               
               ZIss(j,i) = ZIss(j,i) + ZIs(j1,i1) .*(1./(sqrt(2.*pi).*sigma)).* exp(-(1./2).*((sqrt((i-i1).^2+(j-j1).^2)).^2)./(sigma.^2));
               normgauss = normgauss + (1./(sqrt(2.*pi).*sigma)).* exp(-(1./2).*((sqrt((i-i1).^2+(j-j1).^2)).^2)./(sigma.^2)); 
           end
        end  
        
        ZIss(j,i) = ZIss(j,i)./normgauss;
        
    end
end

%H=fspecial('gaussian',zeropad,zeropad);
%ZIs = imfilter(ZIs,H,'replicate');

  vp_resamp(1:ny,1:nx)=ZIss(1+zeropad:ny+zeropad,1+zeropad:nx+zeropad);
  clear ZIss;
end
	imagesc(x,y,vp_resamp);
	hold on;
	%contour(x,y,model1,10,'k-');
    %plot(xrec,yrec,'wo');
	%plot(xshot,yshot,'r*');
	
    %caxis([caxis_value_vp1 caxis_value_vp2]);
    
    colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
       axis([min(x)+FW max(x)-FW min(y) max(y)-FW])
		 axis ij
       
       xlabel('x [m]');
       ylabel('y [m]');
       
       iter_text=['P-Wave Velocity [m/s]'];
       
       title(iter_text);

if(WRITEMODE==1)       
file1=['marmousi_II_scaled_smooth.vp'];
fid1=fopen(file1,'w','ieee-le');
fwrite(fid1,vp_resamp,'float')
fclose(fid1);
end

clear x;
clear y;

% -------------------------------------------------------------------------
% S-Wave Velocity
% -------------------------------------------------------------------------
 
figure;
% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=13601;
NY1=1; NY2=2801; 
IDX=1; IDY=1;
dh=1.25;


% time increment for snapshots:
TSNAPINC=2.0e-2; TSNAP1=0.001;
FW=0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

%load 'seismic.map'
colormap(jet(256));

% load model
 file='Model_vs.bin';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 vs=fread(fid,[ny,nx],'float');
 fclose(fid);
 
 % cut water column
 vs_wo_water = vs(246:2800,:);
 vs_wo_water(1:160,:) = 303.4;
 clear vs;
 
 size(vs_wo_water)
 
 % resample model
 vs_resamp = vs_wo_water(1:2:2475,1:2:13600); 
 
 clear vs_wo_water;
 clear x;
 clear y;

vs_resamp_1 = vs_resamp(1:2:1160,2801:2:4800);
clear vs_resamp;
vs_resamp = vs_resamp_1;
 
 size(vs_resamp) 
 
% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=1000;
NY1=1; NY2=580; 
IDX=1; IDY=1;
dh=5.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

vs_resamp = vp_resamp ./ sqrt(3);

vs_resamp(1:water_depth,:) = 0.0;

caxis_value_vs1 = min(min(vs_resamp)) 
caxis_value_vs2 = max(max(vs_resamp))

if(SMOOTH==1)
% apply median filter
% add boundary to the model to avoid filter artefacts on the inversion grid 
ZIs = ones(ny+2.*zeropad,nx+2.*zeropad);

for j=1:ny
ZIs(j+zeropad,1:zeropad) = vs_resamp(j,1);
ZIs(j+zeropad,nx+zeropad:nx+2.*zeropad) = vs_resamp(j,1);
end

ZIs(1+zeropad:ny+zeropad,1+zeropad:nx+zeropad) = vs_resamp(1:ny,1:nx); 

for i=1:nx+(2.*zeropad)
ZIs(1:zeropad,i) = ZIs(1+zeropad,i);
ZIs(ny+zeropad:ny+2.*zeropad,i) = ZIs(ny+zeropad,i);
end

% apply median filter
  sigma = zeropad;

for j=1:ny+2.*zeropad
    for i=1:nx+2.*zeropad
       ZIss(j,i) = ZIs(j,i);
    end
end
  
% apply gaussian filter
for j=1+zeropad+filter_depth:ny+zeropad
    for i=1+zeropad:nx+zeropad
        
        ZIss(j,i)=0.0;
        normgauss=0.0;
        
        for j1=j-zeropad:j+zeropad
           for i1=i-zeropad:i+zeropad
               
               ZIss(j,i) = ZIss(j,i) + ZIs(j1,i1) .*(1./(sqrt(2.*pi).*sigma)).* exp(-(1./2).*((sqrt((i-i1).^2+(j-j1).^2)).^2)./(sigma.^2));
               normgauss = normgauss + (1./(sqrt(2.*pi).*sigma)).* exp(-(1./2).*((sqrt((i-i1).^2+(j-j1).^2)).^2)./(sigma.^2)); 
           end
        end  
        
        ZIss(j,i) = ZIss(j,i)./normgauss;
        
    end
end

  vs_resamp(1:ny,1:nx)=ZIss(1+zeropad:ny+zeropad,1+zeropad:nx+zeropad);
  clear ZIs;
end

	imagesc(x,y,vs_resamp);
	hold on;
	%contour(x,y,model1,10,'k-');
    %plot(xrec,yrec,'wo');
	%plot(xshot,yshot,'w*');
	
    caxis([caxis_value_vs1 caxis_value_vs2]);
    
    colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
       axis([min(x)+FW max(x)-FW min(y) max(y)-FW])
		 axis ij
       
       xlabel('x [m]');
       ylabel('y [m]');
       
       iter_text=['S-Wave Velocity [m/s]'];
       
       title(iter_text);
       
       hold off;

% vs_resamp(1,:) = 1e-6;
% vs_resamp(2,:) = 1e-6;

if(WRITEMODE==1) 
file1=['marmousi_II_scaled_smooth.vs'];
fid1=fopen(file1,'w','ieee-le');
fwrite(fid1,vs_resamp,'float')
fclose(fid1);
end

clear x;
clear y;
clear vs_resamp;

% -------------------------------------------------------------------------
% Density
% -------------------------------------------------------------------------

% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=13601;
NY1=1; NY2=2801; 
IDX=1; IDY=1;
dh=1.25;

FW=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

%load 'seismic.map'
colormap(jet(256));

% load model
 file='Model_rho.bin';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 rho=fread(fid,[ny,nx],'float');
 fclose(fid);
 
 % cut water column
 rho_wo_water = rho(246:2800,:);
 rho_wo_water(1:160,:) = 1.956;
 
 clear rho;
 
 size(rho_wo_water)
 
 % resample model
 rho_resamp = rho_wo_water(1:2:2475,1:2:13600).*1000.0; 
 
 clear rho_wo_water;
 clear x;
 clear y;
 
rho_resamp_1 = rho_resamp(1:2:1160,2801:2:4800);
clear rho_resamp;
rho_resamp = rho_resamp_1;
 
 size(rho_resamp)

caxis_value_rho1 = min(min(rho_resamp)) 
caxis_value_rho2 = max(max(rho_resamp))

 
% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=1000;
NY1=1; NY2=580; 
IDX=1; IDY=1;
dh=5.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

rho_resamp = 1000.0 .* 0.31 .* (vp_resamp).^(1./4);

rho_resamp(1:water_depth,:) = 1000.0;

caxis_value_rho1 = min(min(rho_resamp)) 
caxis_value_rho2 = max(max(rho_resamp))

if(SMOOTH==1)
% apply median filter
% add boundary to the model to avoid filter artefacts on the inversion grid 
ZIs = ones(ny+2.*zeropad,nx+2.*zeropad);

for j=1:ny
ZIs(j+zeropad,1:zeropad) = rho_resamp(j,1);
ZIs(j+zeropad,nx+zeropad:nx+2.*zeropad) = rho_resamp(j,1);
end

ZIs(1+zeropad:ny+zeropad,1+zeropad:nx+zeropad) = rho_resamp(1:ny,1:nx); 

for i=1:nx+(2.*zeropad)
ZIs(1:zeropad,i) = ZIs(1+zeropad,i);
ZIs(ny+zeropad:ny+2.*zeropad,i) = ZIs(ny+zeropad,i);
end

% apply median filter
  sigma = zeropad;

for j=1:ny+2.*zeropad
    for i=1:nx+2.*zeropad
       ZIss(j,i) = ZIs(j,i);
    end
end
  
  
% apply gaussian filter
for j=1+zeropad+filter_depth:ny+zeropad
    for i=1+zeropad:nx+zeropad
        
        ZIss(j,i)=0.0;
        normgauss=0.0;
        
        for j1=j-zeropad:j+zeropad
           for i1=i-zeropad:i+zeropad
               
               ZIss(j,i) = ZIss(j,i) + ZIs(j1,i1) .*(1./(sqrt(2.*pi).*sigma)).* exp(-(1./2).*((sqrt((i-i1).^2+(j-j1).^2)).^2)./(sigma.^2));
               normgauss = normgauss + (1./(sqrt(2.*pi).*sigma)).* exp(-(1./2).*((sqrt((i-i1).^2+(j-j1).^2)).^2)./(sigma.^2)); 
           end
        end  
        
        ZIss(j,i) = ZIss(j,i)./normgauss;
        
    end
end

  rho_resamp(1:ny,1:nx)=ZIss(1+zeropad:ny+zeropad,1+zeropad:nx+zeropad);
  clear ZIs;
end

figure; 
	imagesc(x,y,rho_resamp);
    
    hold on;
	%contour(x,y,model1,10,'k-');
    %plot(xrec,yrec,'wo');
	%plot(xshot,yshot,'w*');
	
    caxis([caxis_value_rho1 caxis_value_rho2]);
    
    colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
       axis([min(x)+FW max(x)-FW min(y) max(y)-FW])
		 axis ij
       
       xlabel('x [m]');
       ylabel('y [m]');
       
       iter_text=['Density [kg/m^3]'];
       
       title(iter_text);

if(WRITEMODE==1)       
% rho_resamp(1,:) = 1.25;
% rho_resamp(2,:) = 1.25;

file1=['marmousi_II_scaled_smooth.rho'];
fid1=fopen(file1,'w','ieee-le');
fwrite(fid1,rho_resamp,'float')
fclose(fid1);
end




















