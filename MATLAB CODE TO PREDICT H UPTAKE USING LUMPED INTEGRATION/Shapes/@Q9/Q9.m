classdef Q9
    %Q9 9 node, square, quadratic surface element
    
    properties
        ipcount1D % number of integration points per direction
        ipcount % total number of integration points
        
        rectangular % element is strictly rectangular
        
        Nbase   % basis functions; ip, shapefunc
        Gbase   % basis function derivatives; ip, shapefunc, dx/dy
        
        wbase   % integration point weights; ip
        xbase   % ntegration point coordiantes; ip
        ybase   % ntegration point coordiantes; ip
    end
    
    methods
        function obj = Q9(ipcount1D, rect, zeroWeight)
			% Initializations 
            [x1D, w1D] = obj.getIpscheme(ipcount1D, zeroWeight);
            
            if (zeroWeight)
                ipcount1D = ipcount1D+2;
            end
            obj.ipcount1D = ipcount1D;
            obj.ipcount = obj.ipcount1D*obj.ipcount1D;
            obj.rectangular = rect;
            
            k = 0;
            for j=1:obj.ipcount1D
                for i=1:obj.ipcount1D
                    k=k+1;
                    xbase(k) = x1D(i);
                    ybase(k) = x1D(j);
                    wbase(k) = w1D(i)*w1D(j);
                end
            end
            
            for i=1:length(xbase)
                x = xbase(i); y = ybase(i);
                x1 = -1; x2 = 0; x3 = 1;
                Nx = [(x-x2)*(x-x3)/((x1-x2)*(x1-x3)); 
                      (x-x1)*(x-x3)/((x2-x1)*(x2-x3)); 
                      (x-x1)*(x-x2)/((x3-x1)*(x3-x2))];
                Dx = [1*(x-x3)/((x1-x2)*(x1-x3))+(x-x2)*1/((x1-x2)*(x1-x3)); 
                      1*(x-x3)/((x2-x1)*(x2-x3))+(x-x1)*1/((x2-x1)*(x2-x3)); 
                      1*(x-x2)/((x3-x1)*(x3-x2))+(x-x1)*1/((x3-x1)*(x3-x2))];
                Ny = [(y-x2)*(y-x3)/((x1-x2)*(x1-x3)); 
                      (y-x1)*(y-x3)/((x2-x1)*(x2-x3)); 
                      (y-x1)*(y-x2)/((x3-x1)*(x3-x2))];
                Dy = [1*(y-x3)/((x1-x2)*(x1-x3))+(y-x2)*1/((x1-x2)*(x1-x3)); 
                      1*(y-x3)/((x2-x1)*(x2-x3))+(y-x1)*1/((x2-x1)*(x2-x3)); 
                      1*(y-x2)/((x3-x1)*(x3-x2))+(y-x1)*1/((x3-x1)*(x3-x2))];
                
                Nbase(i,:) = kron(Ny, Nx);
                Gbase(i,:,1) = kron(Ny, Dx);
                Gbase(i,:,2) = kron(Dy, Nx);
            end
            
            obj.xbase = xbase;
            obj.ybase = ybase;
            obj.wbase = wbase;
            obj.Nbase = Nbase;
            obj.Gbase = Gbase;
            
        end
        
        function [N, G, w] = getVals(obj, X, Y)
			% Get shape function values, gradients, and accompanying
			% integration weights, providing the nodal coordinate vectors X
			% and Y

            N = obj.Nbase;
            
            if obj.rectangular 
            	G(:,:,1) = obj.Gbase(:,:,1)*2/(X(3)-X(1));
            	G(:,:,2) = obj.Gbase(:,:,2)*2/(Y(7)-Y(1));
                w = obj.wbase * (X(3)-X(1))/2 * (Y(7)-Y(1))/2;
            else 
                dXdXi = obj.Gbase(:,:,1)*X; dXdEta = obj.Gbase(:,:,2)*X;
                dYdXi = obj.Gbase(:,:,1)*Y; dYdEta = obj.Gbase(:,:,2)*Y;

                J(:,1,1) = dXdXi; J(:,1,2) = dXdEta;
                J(:,2,1) = dYdXi; J(:,2,2) = dYdEta;

                G = 0.0*obj.Gbase;
                for i=1:size(obj.Nbase, 1)
                    Jinv = inv(squeeze(J(i,:,:)));
                    for j=1:size(obj.Nbase, 2)
                        G(i,j,:) = Jinv*squeeze(obj.Gbase(i,j,:));
                    end

                    w = obj.wbase(i)*det(squeeze(Jinv));
                end
            
            end

        end
        
        function xy = getIPGlobal(obj, X,Y)
			% Get global coordinates for integration points
            
            if obj.rectangular 
                xy(1,:) = X(1) + (0.5*(obj.xbase+1))*(X(3)-X(1));
                xy(2,:) = Y(1) + (0.5*(obj.ybase+1))*(Y(7)-Y(1));
            else
                %% not implemented
            end

        end
 
    end
    
    
    methods (Access = private)
        function [x1D, w1D] = getIpscheme(obj, ipcount1D, zeroWeight)
			% Gauss integration scheme
            if (ipcount1D == 1)
                x1D = 0;
                w1D = 2;
            elseif (ipcount1D == 2)
                x1D = [-1/sqrt(3); 1/sqrt(3)];
                w1D = [1; 1];                
            elseif (ipcount1D == 3)
                x1D = [-sqrt(3/5); 0; sqrt(3/5)];
                w1D = [5/9; 8/9; 5/9];                     
            elseif (ipcount1D ==4)
                x1D = [-sqrt(3/7+2/7*sqrt(6/5)); -sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7+2/7*sqrt(6/5))];
                w1D = [(18-sqrt(30))/36; (18+sqrt(30))/36; (18+sqrt(30))/36; (18-sqrt(30))/36];                     
            elseif (ipcount1D == 5)
                x1D = [-1/3*sqrt(5+2*sqrt(10/7)); -1/3*sqrt(5-2*sqrt(10/7)); 0; 1/3*sqrt(5-2*sqrt(10/7)); 1/3*sqrt(5+2*sqrt(10/7))];
                w1D = [(322-13*sqrt(70))/900; (322+13*sqrt(70))/900; 128/225; (322+13*sqrt(70))/900; (322-13*sqrt(70))/900];                     
            else
                error("Higer order ip schemes not implemented in Shapes.Q9");
            end
            
            if (zeroWeight)
               x1D = [-1; x1D; 1];
               w1D = [0; w1D; 0];
            end
            
        end
    end
end

