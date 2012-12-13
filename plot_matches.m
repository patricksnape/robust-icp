function plot_matches(A, B)

    hndl = figure; % creates a plotting window and stores the handle in hndl
    set(hndl,'Renderer','OpenGL');
    hold on;
    scatter(A, 'r.');  
    scatter(B, 'yo');

    N = size(A, 1);

    for i = 1:N
        line([A(i, 1), B(i, 1)], ...
             [A(i, 2), B(i, 2)], ...
             [A(i, 3), B(i, 3)]);
    end
    hold off;
end