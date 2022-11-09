function edgecoords = SPIE_JMI_hand_findedgecoords(I)
    % determine the correct coordinates of the ONLY segment in the data
    % (there shouldn't be any holes in the binarized data)
    I = double(I);
    [B,L] = bwboundaries(imbinarize(I));
    edgecoords = B{1}; edgecoords = [edgecoords(:,2) edgecoords(:,1)];
end