function plotJvsTheta(obj)

if isempty(obj.theAx)
    figure;
    obj.theAx = axes;
else
    axes(obj.theAx);
end

% surf
if size(obj.theta, 2) == 2
    hold on
    tri = delaunay(obj.theta(:,1), obj.theta(:,2));
    theSurf = trisurf(tri, obj.theta(:,1), obj.theta(:,2), ...
        obj.objectives);
    theSurf.FaceAlpha = 0.75;
    theSurf.LineStyle = ':';
    theSurf.EdgeAlpha = 0.75;
    view(-60, 30)
    xlabel('\theta_1')
    ylabel('\theta_2')
    zlabel('objective')
    title(['Objective ' obj.name])
else
    % parallelcoord
    Y = sortrows([obj.objectives obj.theta],1);
    N = size(Y,1);
    c = jet(N);
    % standardize
    Y = (Y - repmat(min(Y), N, 1)) ./ repmat(max(Y) - min(Y), N, 1);
    labels = {'Objective'};
    for th = 1:size(obj.theta,2)
        labels{end+1} = ['Theta #' num2str(th)];
    end
    hold on
    for i = 1:N
        parallelcoords(Y(i,:), 'Color', c(i,:), 'Labels', labels);
    end
    ylabel('Normalized values')
    title(['Objective: ' obj.name])
end