function dt = pointclouddt(cloud, N)

%%% N is the size of the volume NxNxN

B = round(cloud);
M = zeros(N,N,N);
M(sub2ind(size(M), B(:,1), B(:,2), B(:,3))) = 1;
dt = bwdistX(M);

end