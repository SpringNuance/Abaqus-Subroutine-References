function [Nodes, Elementgroups, Nodegroups, Area, rect] = Geometry_Generator(obj, props)
	%GEOMETRY_GENERATOR Defines the geometry and mesh used within the
	%simulations. Inputs required:
	%   mesh_in.Lx    = 10e-3;	%Domain length [m]
    %	mesh_in.Ly    = 10e-3;	%Domain Height [m]
    %	mesh_in.ipcount1D = 3;	%numer of integration points per dimension
    %	mesh_in.zeroWeight = true; %Include zero-weight integration points at element boundaries (for post-processing)
	
	if true %if true, generate anew, if false then load from file
    	Lx = props.Lx; 
    	Ly = props.Ly;  
	
		%metal
		R1 = [3,4,0,Lx,Lx,0,0,0,Ly,Ly]';
		C1 = [1,5e-3-0.2e-3,5e-3,0.2e-3]';
		C1 = [C1;zeros(length(R1) - length(C1),1)];
		R2 = [3,4,0,5e-3-0.2e-3,5e-3-0.2e-3,0,Ly/2-0.2e-3,Ly/2-0.2e-3,Ly/2+0.2e-3,Ly/2+0.2e-3]';
		R3 =  [3,4,-Lx,0,0,-Lx,0,0,Ly,Ly]';
		gm = [R1,C1,R2,R3];
		sf = '(R1-C1-R2)+(R3+C1+R2)';
	
		ns = char('R1','C1','R2','R3');
		ns = ns';
		[shp, shpb] = decsg(gm,sf,ns);
		[shp, shpb] = csgdel(shp,shpb,[3, 9, 13, 16]);
	
	% 	subplot(2,1,1)
	% 	pdegplot(shp,'EdgeLabels','on','FaceLabels','on')
	
	
		geo = createpde(1);
		geometryFromEdges(geo,shp);
		generateMesh(geo,'Hmax',1e-3,'Hgrad',1.2,'Hedge',{[2,3,7,8,11,12], 0.1e-3});
	
	% 	subplot(2,1,2)
	% 	pdeplot(geo,'NodeLabels','off','ElementLabels','off')
	
		Nodes = geo.Mesh.Nodes';
	
		% interior elements
		Elementgroups{1}.name = "Metal";
		Elementgroups{1}.type = "T6";
		Elementgroups{1}.Elems = geo.Mesh.Elements(:,findElements(geo.Mesh,'region','Face',1))';
	
		Elementgroups{2}.name = "Electrolyte";
		Elementgroups{2}.type = "T6";
		Elementgroups{2}.Elems = geo.Mesh.Elements(:,findElements(geo.Mesh,'region','Face',2))';
	
		%% exterior boundary elements
		Elementgroups{3}.name = "E_Left";
		Elementgroups{3}.type = "L3B";
		
		N = findNodes(geo.Mesh,'region','Edge',[4]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,2));
	
		cntr = 0;
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{3}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
	
	
	
		Elementgroups{4}.name = "M_Right";
		Elementgroups{4}.type = "L3B";
		
		N = findNodes(geo.Mesh,'region','Edge',[1]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,2));
	
		cntr = 0;
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{4}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
	
	
		Elementgroups{5}.name = "E_Bottom";
		Elementgroups{5}.type = "L3B";
		
		N = findNodes(geo.Mesh,'region','Edge',[5]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,1));
	
		cntr = 0;
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{5}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
	
	
		Elementgroups{6}.name = "E_Top";
		Elementgroups{6}.type = "L3B";
		
		N = findNodes(geo.Mesh,'region','Edge',[9]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,1));
	
		cntr = 0;
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{6}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
	
	
	
		Elementgroups{7}.name = "M_Bottom";
		Elementgroups{7}.type = "L3B";
		
		N = findNodes(geo.Mesh,'region','Edge',[6]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,1));
	
		cntr = 0;
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{7}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
	
	
		Elementgroups{8}.name = "M_Top";
		Elementgroups{8}.type = "L3B";
		
		N = findNodes(geo.Mesh,'region','Edge',[10]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,1));
	
		cntr = 0;
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{8}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
	
	
		%% metal-electrolyte boundary
		Elementgroups{9}.name = "Interface";
		Elementgroups{9}.type = "L3B";
	
		[p,e,t] = meshToPet(geo.Mesh);
	
		N = findNodes(geo.Mesh,'region','Edge',[7, 2, 11, 12, 3, 8]);
		xy = Nodes(N,:);
		[~,ind] = min(xy(:,2));
	
		stop = false;
		cnt = 0;
		N0 = N(ind);
		while stop==false
			cnt=cnt+1;
			nds = N0;
	
			els = findElements(geo.Mesh,'attached',N0);
			ns = [];
			for i=1:length(els)
				ns = [ns;geo.Mesh.Elements(:,els(i))];
			end
			ns = unique(ns);
			bc = [];
			for i=1:length(ns)
				z = find(ns(i)==N);
				if (~isempty(z))
					bc = [bc ns(i)];
				end
			end
			if length(bc)>3
				xy0 = Nodes(N0,:);
				xy = Nodes(bc,:);
				for k=1:length(bc)
					dst(k) = (xy(k,1)-xy0(1)).^2+(xy(k,2)-xy0(2)).^2;
				end
				[~,indc] = mink(dst,3);
				bc = bc(indc);
			end
			if (length(bc)==1)
				stop = true; %error here
			else
				bc = bc(bc~=nds(1));
				nds(2) = max(bc);
				nds(3) = min(bc);
		
				N=N(N~=nds(1));
				N=N(N~=nds(2));
				N0 = nds(3);
		
				Elementgroups{9}.Elems(cnt,:) = nds;
				if (length(N)==1)
					stop=true;
				end
			end
			
		end
	
    	for g=1:length(Elementgroups)
        	Nodegroups{g}.name = Elementgroups{g}.name;
        	Nodegroups{g}.Nodes = unique(reshape(Elementgroups{g}.Elems,[],1));
    	end
    	
		Area = zeros(9,1);
		rect = false;

		save('mesh.mat','Nodes', 'Elementgroups', 'Nodegroups', 'Area', 'rect')
	else
		load('mesh.mat','Nodes', 'Elementgroups', 'Nodegroups', 'Area', 'rect')
	end
end

