% ===================================================
% picking first arrivals for shallow seismic datasets
% contact: m.schaefer (at) kit.edu - this is just
% to get an idea how it could be done
% note: this matlab script requires 'su2matlab'
% ===================================================

close all
clear all

% short READE ME
% 1) please choose your input file (SU-format) -> filename_in
% 2) choose geophon distant -> geop_dist (m)
% 3) it is possible to apply a time window in which the trigger
%    is applied (tw_start, tw_end)
% 4) vector with triggering points (samples) will be saved -> filename_out
%
% !please note that only the first picks of every trace are plotted and saved!
%
% the most important parameters for the sta/lta picker are:
% STA      Short-term averaging window length (s)
% LTA      Long-term averaging window length (s)
% TR       Value of STA/LTA ratio that triggers


% input
filename_in='../su/measured_data/rheinstetten_24042012_y.su.shot1';
data=su2matlab(filename_in);
fileparameters = dir(filename_in);
no_of_traces=fileparameters.bytes/(240+data(1).ns*4);

% geophon distant (m)
geop_dist=1.5;
h=(1:no_of_traces)*geop_dist;

% end of time window for trigger (samples)
tw_start=1;
tw_end=2000;

% output
filename_out='picks_1.dat';

maxx=zeros(no_of_traces);
picked_times{:}=zeros(no_of_traces);

for i=1:no_of_traces
    maxx(i)=max(abs(data(i).trace));
    data(i).trace=data(i).trace/maxx(i)+h(i);
    
    % sta/lta Picker
    sig=data(i).trace(tw_start:tw_end);
    DT=data(i).dt*10^(-6);
    BE=[0 0.7];
    STA=0.007;
    LTA=0.08;
    TR=2;
    DTR=1;
    PEM=0;
    PET=0;
    PNL=0;
    ATL=0;
    
    [trigt,stav,ltav,ratio,tim1,tim2,tim3,trigs,dtrigs]=stalta(sig,DT,BE,STA,LTA,TR,DTR,PEM,PET,PNL,ATL);
    
    % picked times
    picked_times{i}=trigs*data(1).dt*10^(-6);
    
    % sample/time vector
    x=(1:data(1).ns)*data(1).dt*10^(-6);
    
    figure(1)
    hold on
    plot(x,data(i).trace)
    plot(picked_times{i}(1,1),h(i),'+r','MarkerSize',12)
    %xlim([0 0.7])
    ylim([0 73])
    xlabel('time (s)')
    ylabel('offset (m)')
    title('Normalized seismograms with first triggering points')
    set(get(gca,'title'),'FontSize',12);
    set(get(gca,'title'),'FontWeight','bold');
    set(get(gca,'Ylabel'),'FontSize',12);
    set(get(gca,'Ylabel'),'FontWeight','bold');
    set(get(gca,'Xlabel'),'FontSize',12);
    set(get(gca,'Xlabel'),'FontWeight','bold');
    set(gca,'FontSize',12);
    set(gca,'FontWeight','bold');
    set(gca,'Linewidth',1.0);
    set(gca,'Box','on');
end

fid = fopen(filename_out,'w');
for i=1:no_of_traces
fprintf(fid,'%d ',picked_times{i}(1,1));
fprintf(fid,'\n');
end
fclose(fid);
