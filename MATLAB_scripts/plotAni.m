function plotAni(spatialX, spatialY)
    global fAni
    try
        close(fAni)
    end
    fAni = figure;

    for i = 1 : size(spatialX,1)
        if ~ishandle(fAni)
            break
        end
        plot(spatialX(i,:),spatialY(i,:),'Color','black')
        axis equal
        xlim([min(spatialX(:)),max(spatialX(:))])
        ylim([min(spatialY(:)),max(spatialY(:))])
        drawnow
        pause(0.05)
    end
end