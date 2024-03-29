function PlotIP(obj, varName, plotloc)
	%PLOTIP Plots integration-point level data

    for g=1:length(obj.mesh.Elementgroups)
        if (obj.mesh.Elementgroups{g}.name == plotloc)
            for e=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                IPvar = obj.Request_Info(varName, e, "Interior");
                xy = obj.mesh.getIPCoords(g, e);
                
                x(e, :) = xy(1,:);
                y(e, :) = xy(2,:);
                z(e, :) = IPvar;
            end
            
            xlims = [min(min(x)), max(max(x))];
            ylims = [min(min(y)), max(max(y))];
            
            [xi,yi] = meshgrid(linspace(xlims(1), xlims(2), 1000), linspace(ylims(1), ylims(2), 1000));
            zi = griddata(x,y,z,xi,yi,'linear');
            s = pcolor(xi,yi,zi);
            s.EdgeColor = 'none';
            hold on
            colorbar
        end
    end
end


