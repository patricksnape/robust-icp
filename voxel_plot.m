function voxel_plot(mat)

%VOXELPLOT function to draw a 3-D binary image
%
%Usage
%   voxelPlot(mat);
%
%   Draws a 3D binary image
%   Relies on Suresh Joel's voxel function
%
%   mat is an NxNxN binary matrix of type single
%

%   Patrick Snape 04 Dec 2012

ind = find(mat);
for i = 1:size(ind)
    [x, y, z] = ind2sub(size(mat), ind(i));
    voxel([x, y, z]);
end

end