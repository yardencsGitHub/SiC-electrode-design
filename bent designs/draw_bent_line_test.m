function [S,X,Y,Ytag,K] = draw_bent_line_test(x0,y0,y_tag0,k0,dk_ds,ds,target_s)
    figure;
    axes; hold on;
    chan_centers = y0;
    y=y0; x=x0; y_tag=y_tag0; k=k0; s=0;
    S=[0]; X=[x]; Y=[y]; Ytag=[y_tag]; K=[k];
    while (s < target_s)
        r=abs(1/k); phi=abs(ds/r);
        beta = abs(atan(y_tag));
        y_end = y + 2*r*sin(phi/2)*sin(beta)*sign(chan_centers-0.5);
        x_end = x + 2*r*sin(phi/2)*cos(beta);
        x_center = x + r*sin(phi/2-beta); 
        y_center = y + r*cos(phi/2-beta)*sign(chan_centers-0.5);
        t0=asin((y-y_center)/r);
        t1=t0+phi*sign(chan_centers-0.5);
        t=t0:abs(t1-t0)/10*sign(chan_centers-0.5):t1;
        plot(x_center+r*cos(t),y_center+r*sin(t),'k');
        s = s + ds;
        k = k + dk_ds;
        y_tag = tan(beta + phi); x = x_end; y = y_end;
        S = [S;s]; Ytag = [Ytag;y_tag]; X = [X;x_end]; Y = [Y;y_end]; K = [K;k];
    end
end