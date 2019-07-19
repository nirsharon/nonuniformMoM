% script name: "test_generateSOS_B_forProvable"

clear;
generateSOS_B_forProvable;

S       = load('SO3_fifteen.mat');
SO_grid = S.SO3;
P_B = numel(B);

% positiveity check
vals = zeros(size(SO_grid,3),1);
for i=1:size(SO_grid,3)
    eul     = rotm2eul(SO_grid(:,:,i));
    vals(i) = wignerD_expansion2(B,eul);
end
v = min(real(vals));
fprintf('The minimum value over the grid is %d \n', v)
