function xi = SPIEJMI_getXiStars(arc,n)
    % find the special points in an arc based on the number of subintervals
    % specified
    arclength = size(arc,1);
    deltaX = arclength/n;
    x_ind = round(1:deltaX:arclength);
    xi = arc(x_ind,:);
end  