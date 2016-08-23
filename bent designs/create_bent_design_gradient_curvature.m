function create_bent_design_gradient_curvature(total_length,bent_start,pad_size,chan_num,chan_width,chan_space,dk_ds_max)
ds=1e-4;
% This script creates a curved electrode assembly design with the following
% geometry and parameters:
% The straight shank length (in mm) is set by the input parameter
% 'total_length' and used to normalize all other measurements. The
% resulting device is scaled to 1 = total_length.
% The design is made of shanks that are evenly apaced at their base, start with a straight segment,
% and spread in curvatures that increase in a constant rate. 
% The device parameters:
%   chan_num    - Number of shanks.
%   bent_start  - The length (in mm) of the straight stem of the shanks.
%   pad_size    - The length (in mm) of the electrode pads at the tips.
%   chan_width  - The width of each shank (in mm).
%   chan_space  - The spacing between shanks (in mm) at the base of the device.
%   dk_ds_max   - The maximal curving gradient of the device. (in 1/radians/length)


% total_length=5;
% edge_dim=1;
% bent_start=0.1;
% pad_size=0.1;
% chan_num=16;
% chan_width=0.014;
% chan_space=0.025;
% phi_max=5/180*pi;

    pad_width_ratio = 0.002 / total_length;
    edge_point_ratio = 0.02 / total_length;
    width_ratio = chan_width / total_length;
    straight_ratio = bent_start / total_length;
    pad_ratio = pad_size / total_length;
    %edge_ratio = edge_dim / total_length;
    spacing_ratio = chan_space / total_length;
    chan_centers = flipud(0.5 + spacing_ratio*[-(chan_num-1)/2:(chan_num-1)/2]');
    %end_centers = flipud(0.5 + edge_ratio*[-(chan_num-1)/2:(chan_num-1)/2]'/(chan_num-1));
    ddkds = 2*dk_ds_max / (chan_num - 1);
    delta_k = 0*ddkds / spacing_ratio / 20; % curvature gradient's gradient
    dkds = abs(-dk_ds_max + [0:(chan_num-1)]'*ddkds);

    h = figure; axis; hold on;
    Ss=[]; Xs=[]; Ys=[]; Ytags=[]; Ks=[];
    for chan = 1:chan_num
        [S1,X1,Y1,Ytag1,K1] = draw_bent_line(chan_centers(chan),straight_ratio,chan_centers(chan)+width_ratio/2,0,1e-6,dkds(chan)+delta_k*width_ratio/2*sign(chan_centers(chan)-0.5),ds,1,[0 width_ratio*(chan_centers(chan)>0.5)],1);  %0*
        line([0 min(X1)],[chan_centers(chan)-width_ratio/2 chan_centers(chan)-width_ratio/2],'Color','k');
        x1 = X1(end); y1 = Y1(end);
        [S,X,Y,Ytag,K] = draw_bent_line(chan_centers(chan),straight_ratio,chan_centers(chan)-width_ratio/2,0,1e-6,dkds(chan)-delta_k*width_ratio/2*sign(chan_centers(chan)-0.5),ds,1,[0 width_ratio*(chan_centers(chan)<=0.5)],1); %-0*
        line([0 min(X)],[chan_centers(chan)+width_ratio/2 chan_centers(chan)+width_ratio/2],'Color','k');
        tipx = (X(end)+x1)/2 + cos(atan(Ytag(end)))*edge_point_ratio; 
        tipy = (Y(end)+y1)/2 + sin(atan(Ytag(end)))*edge_point_ratio*sign(chan_centers(chan)-0.5);
        line([x1 tipx],[y1 tipy],'Color','k');
        line([X(end) tipx],[Y(end) tipy],'Color','k');
        [S,X,Y,Ytag,K] = draw_bent_line(chan_centers(chan),straight_ratio,chan_centers(chan),0,1e-6,dkds(chan),ds,1,[0 0],0);
        pad_start_ind = max(find(S < (1-pad_ratio-width_ratio/2*abs(Ytag(end)))));
        x0 = X(pad_start_ind)-pad_width_ratio/2*sign(chan_centers(chan)-0.5)*abs(sin(atan(Ytag(pad_start_ind))));
        y0 = Y(pad_start_ind)+pad_width_ratio/2*abs(cos(atan(Ytag(pad_start_ind))));
        %(X(pad_start_ind)+X1(pad_start_ind))/2;
        %y0 = (Y(pad_start_ind)+Y1(pad_start_ind))/2;
        k0 = K(pad_start_ind); ytag0 = Ytag(pad_start_ind);
        [pS,pX,pY,pYtag,pK]=draw_bent_line(chan_centers(chan),x0,y0,ytag0,k0,dkds(chan),ds,pad_ratio,pad_width_ratio/2*(chan_centers(chan)>0.5)*[0 1],1);
        x0 = X(pad_start_ind)+pad_width_ratio/2*sign(chan_centers(chan)-0.5)*abs(sin(atan(Ytag(pad_start_ind))));
        y0 = Y(pad_start_ind)-pad_width_ratio/2*abs(cos(atan(Ytag(pad_start_ind))));
        [S,X,Y,Ytag,K]=draw_bent_line(chan_centers(chan),x0,y0,ytag0,k0,dkds(chan),ds,pad_ratio,pad_width_ratio/2*(chan_centers(chan)<=0.5)*[0 1],1);
        line([pX(1) X(1)],[pY(1) Y(1)],'Color','k');
        line([pX(end) X(end)],[pY(end) Y(end)],'Color','k');
    end
    function [S,X,Y,Ytag,K] = draw_bent_line(chan_center,x0,y0,y_tag0,k0,dk_ds,ds,target_s,Srange,to_draw) %updown
        y=y0; x=x0; y_tag=y_tag0; k=k0; s=0;
        S=[0]; X=[x]; Y=[y]; Ytag=[y_tag]; K=[k];
        while  (s <= target_s) %(target_s-updown*abs(y_tag)) (s <= target_s-Srange(2)*y_tag) %
            r=abs(1/k); phi=abs(ds/r);
            beta = abs(atan(y_tag));
            y_end = y + 2*r*sin(phi/2)*sin(beta)*sign(chan_center-0.5);
            x_end = x + 2*r*sin(phi/2)*cos(beta);
            x_center = x + r*sin(phi/2-beta); 
            y_center = y + r*cos(phi/2-beta).*sign(chan_center-0.5);
            t0=asin((y-y_center)/r);
            t1=t0+phi*sign(chan_center-0.5);
            t=t0:abs(t1-t0)/10*sign(chan_center-0.5):t1; %10
            %X = [X;x_end]; Y = [Y;y_end];
            X = [X;x_center+r*cos(t')]; Y = [Y;y_center+r*sin(t')];
            S = [S;s+[ds/numel(t):ds/numel(t):ds]'];
            Ytag = [Ytag;tan(beta + [phi/numel(t):phi/numel(t):phi]')];
            K = [K;k+dk_ds*[ds/numel(t):ds/numel(t):ds]'];
            %S = [S;s]; Ytag = [Ytag;y_tag];  K = [K;k];
            %plot(x_center+r*cos(t),y_center+r*sin(t),'k');
            s = s + ds;
            k = k + dk_ds*ds;
            y_tag = tan(beta + phi); x = x_end; y = y_end;
         
%             if (s > (target_s-updown*abs(y_tag)))
%                 break;
%             end
        end
        %tmp = -updown*abs(y_tag0);
        ind = find((S >= Srange(1)*y_tag0) & (S <= (target_s-Srange(2)*y_tag)));
        S=S(ind); X=X(ind); Y=Y(ind); Ytag=Ytag(ind); K=K(ind);
        if (to_draw)
            plot(X,Y,'k');
        end
    end

end