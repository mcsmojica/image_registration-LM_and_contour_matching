function LM_fcnval = SPIEJMI_computeLMfcnval_2d(T,R,alpha)
    % computes squared distance between corresponding points, and sums over
    % all points/landmarks
    % computes squared distance between corresponding points, and sums over
    % all points/landmarks
    if alpha == 0
        LM_fcnval = 0;
    else
        LM_fcnval = alpha*sum((T(:,1) - R(:,1)).^2 +(T(:,2) - R(:,2)).^2);
    end
end