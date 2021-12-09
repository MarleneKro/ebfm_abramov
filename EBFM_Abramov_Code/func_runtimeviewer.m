function  func_runtimeviewer(A,io,t,grid,insol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to plot output while running
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (io.runtimeview && mod(t,io.runtimeview_freq)==0)
    runview = struct;
    
    runview.values{1} = insol.shade;
    runview.values{1} = A.SWin;
    runview.title{1} = 'SWin';
    
    if t==1
        figure; 
        hold on;
    end
    temp = nan(grid.Lx,grid.Ly);
    for n=1:length(grid.xind)
        temp(grid.xind(n),grid.yind(n)) = runview.values{1}(n);
    end
    pcolor(grid.x,grid.y,temp);shading flat; axis equal;
    axis tight; h = colorbar;set(get(h,'ylabel'),'String',runview.title{1});
    caxis([0 450]); 
    drawnow
end

end

