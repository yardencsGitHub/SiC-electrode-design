function create_bent_design(total_length,edge_dim,bent_start,pad_size,chan_num,chan_width,chan_space,phi_max)
% This script creates a curved electrode assembly design with the following
% geometry and parameters

% total_length=5;
% edge_dim=1;
% bent_start=0.1;
% pad_size=0.1;
% chan_num=16;
% chan_width=0.014;
% chan_space=0.025;
% phi_max=5/180*pi;

pad_width_ratio = 0.002 / total_length;
edge_point_ratio = 0.01 / total_length;
width_ratio = chan_width / total_length;
straight_ratio = bent_start / total_length;
pad_ratio = pad_size / total_length;
edge_ratio = edge_dim / total_length;
spacing_ratio = chan_space / total_length;
chan_centers = flipud(0.5 + spacing_ratio*[-(chan_num-1)/2:(chan_num-1)/2]');
end_centers = flipud(0.5 + edge_ratio*[-(chan_num-1)/2:(chan_num-1)/2]'/(chan_num-1));
dphi = 2*phi_max / (chan_num - 1);
phis = abs(-phi_max + [0:(chan_num-1)]'*dphi);
radii = (1-straight_ratio)./abs(phis);
pad_phis = pad_ratio./radii;
edge_phis = edge_point_ratio./radii;
x_end = straight_ratio+sqrt((2*radii.*sin(phis/2)).^2-(end_centers-chan_centers).^2);
betas = abs(atan((end_centers-chan_centers)./(x_end-straight_ratio)));
xs = straight_ratio + radii.*sin(phis/2-betas); 
ys = chan_centers + radii.*cos(phis/2-betas).*sign(chan_centers-0.5);

h = figure; axis; hold on;

for chan = 1:chan_num
    %plot(straight_ratio,chan_centers(chan),'o');
    %plot(x_end(chan),end_centers(chan),'o');
    t0=asin((chan_centers(chan)-ys(chan))/radii(chan));
    t1=t0+phis(chan)*sign(chan_centers(chan)-0.5);
    t=t0:abs(t1-t0)/1000*sign(chan_centers(chan)-0.5):t1;
    plot(xs(chan)+(radii(chan)-width_ratio/2)*cos(t),ys(chan)+(radii(chan)-width_ratio/2)*sin(t),'k');
    plot(xs(chan)+(radii(chan)+width_ratio/2)*cos(t),ys(chan)+(radii(chan)+width_ratio/2)*sin(t),'k');
    
    start_x = min(xs(chan)+(radii(chan)-width_ratio/2)*cos(t0),xs(chan)+(radii(chan)-width_ratio/2)*cos(t1));
    line([0 start_x],[chan_centers(chan)-width_ratio/2 chan_centers(chan)-width_ratio/2],'Color','k');
    start_x = min(xs(chan)+(radii(chan)+width_ratio/2)*cos(t0),xs(chan)+(radii(chan)+width_ratio/2)*cos(t1));
    line([0 start_x],[chan_centers(chan)+width_ratio/2 chan_centers(chan)+width_ratio/2],'Color','k');
    
    
    tipx = max(xs(chan)+(radii(chan))*[cos(t+edge_phis(chan)) cos(t-edge_phis(chan))]);
    tipy = max(ys(chan)+(radii(chan))*[sin(t+edge_phis(chan)) sin(t-edge_phis(chan))])*(chan_centers(chan)>0.5)+...
        min(ys(chan)+(radii(chan))*[sin(t+edge_phis(chan)) sin(t-edge_phis(chan))])*(chan_centers(chan)<=0.5);
    line([max(xs(chan)+(radii(chan)-width_ratio/2)*cos(t)) tipx], ...
        [max(ys(chan)+(radii(chan)-width_ratio/2)*sin(t))*(chan_centers(chan)>0.5)+...
        min(ys(chan)+(radii(chan)-width_ratio/2)*sin(t))*(chan_centers(chan)<=0.5) tipy],'Color','k');
    line([max(xs(chan)+(radii(chan)+width_ratio/2)*cos(t)) tipx], ...
        [max(ys(chan)+(radii(chan)+width_ratio/2)*sin(t))*(chan_centers(chan)>0.5)+...
        min(ys(chan)+(radii(chan)+width_ratio/2)*sin(t))*(chan_centers(chan)<=0.5) tipy],'Color','k');
    %t1=t0+(phis(chan)+edge_phis(chan))*sign(chan_centers(chan)-0.5);
    %t=t0:pi/10000000*sign(chan_centers(chan)-0.5):t1;
    %plot(xs(chan)+(radii(chan))*cos(t),ys(chan)+(radii(chan))*sin(t),'--k');
    %line([xs(chan)+(radii(chan)+width_ratio/2*[-1])*max(cos(t)) xs(chan)+(radii(chan))*max(cos(t))],ys(chan)+(radii(chan)+width_ratio/2*[-1])*max(sin(t)),'Color','k');
    
    t0_pad = t1-pad_phis(chan)*sign(chan_centers(chan)-0.5);
    t=t0_pad:abs(t1-t0)/1000*sign(chan_centers(chan)-0.5):t1;
    plot(xs(chan)+(radii(chan)-pad_width_ratio/2)*cos(t),ys(chan)+(radii(chan)-pad_width_ratio/2)*sin(t),'k');
    plot(xs(chan)+(radii(chan)+pad_width_ratio/2)*cos(t),ys(chan)+(radii(chan)+pad_width_ratio/2)*sin(t),'k');
    line(xs(chan)+(radii(chan)+pad_width_ratio/2*[1 -1])*max(cos(t)),ys(chan)+(radii(chan)+pad_width_ratio/2*[1 -1])*...
        (max(sin(t))*(chan_centers(chan)>0.5)+min(sin(t))*(chan_centers(chan)<=0.5)),'Color','k');
    line(xs(chan)+(radii(chan)+pad_width_ratio/2*[1 -1])*min(cos(t)),ys(chan)+(radii(chan)+pad_width_ratio/2*[1 -1])*...
        (min(sin(t))*(chan_centers(chan)>0.5)+max(sin(t))*(chan_centers(chan)<=0.5)),'Color','k');
    
    %line(xs(chan)+(radii(chan)+pad_width_ratio/2*[-1 1])*min(cos(t)),ys(chan)+(radii(chan)+pad_width_ratio/2*[-1 1])*min(sin(t)),'Color','k');
end

set(gca,'XTick',[]);
set(gca,'YTick',[]);
%xlim([0 1+edge_point_ratio]); ylim([0 1]);
axis equal;
%axis off;
    

     
    
    
