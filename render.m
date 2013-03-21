function render(vertices, triangles, normals)

options = [];
if ~exist('normals', 'var')
    options.normals = [];
else
    options.normal = normals;
end
plot_mesh(vertices, triangles, options);