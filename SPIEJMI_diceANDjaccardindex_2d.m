function [d_ind,j_ind] = SPIEJMI_diceANDjaccardindex_2d(A,B)
    % mia mojica
    % A, B are logical arrays
    % 0 = black pixels, 1 = white
    
    A = double(A); B = double(B);
    card_A = sum(A(:)==1); 
    card_B = sum(B(:)==1);
    
    % pixel is in intersection if the value of that pixel in the sum is 0
    % pixel is in union if the value in the sum is 0 or 1 (- not necessary)
    sum_AB = double(A) + double(B);
    card_AintB = sum(sum_AB(:)==2);
    
    d_ind = 2*card_AintB/(card_A + card_B);
    j_ind = card_AintB/(card_A + card_B - card_AintB);
    
    % number of black voxels multiplied by voxelsize count towards the volume
    pixelarea = 1;
    area_A = sum(A(:)==0)*pixelarea;
    area_B = sum(B(:)==0)*pixelarea;
end