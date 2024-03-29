function Solve(obj)
	%SOLVE solves a single time increment for the nonlinear system of
	%equations through a Newton-Raphson procedure

    stop = false;
    it = 0;

	stepnum = size(obj.convergence_log, 1)+1;
 
	%assemble the stiffness matrix, force vector, and apply the constraints
    obj.physics.Assemble();
    obj.physics.Constrain();

    recalc_pre=true;
    En_err0 = -1;
    curr_max_it = obj.maxIt;
    while(stop == false) %loop untill convergence criterion are fulfilled, or max is exceeded
        fprintf("    Solving it:" + string(it) + "      ");
        tsolve = tic;
        
		%matrix preconditioning 
        recalc_pre = true; 
        if (recalc_pre)
            [P,R,C] = equilibrate(obj.physics.K);
            recalc_pre = false;
        end

		%solve linear system
        if true
            d = -R*P*obj.physics.fint;
            B = R*P*obj.physics.K*C;
			if true %if true, using direct solver, else using iterative gmres
				dy = B\d;
			else
				[L,U] = ilu(B,struct('type','nofill'));
				dy = gmres(B,d,[],1e-4,500,L,U);
			end
            dx = C*dy;
        else
            dx = -obj.physics.K\obj.physics.fint;
        end
        tsolve = toc(tsolve);
        fprintf("        (Solver time:"+string(tsolve)+")\n");

		%perform linear line-search, otherwise just add increment to state
		%vector
        if (obj.linesearch && it>0)
            e0 = obj.physics.fint'*dx;
            obj.physics.Update(dx);
            
            obj.physics.Assemble();
            obj.physics.Constrain();
            
            e1 = obj.physics.fint'*dx;
            factor = -e0/(e1-e0);
            factor = max(obj.linesearchLims(1), min(obj.linesearchLims(2), factor));
            obj.physics.Update(-(1-factor)*dx);
            fprintf("    Linesearch: " + string(e0) + " -> " + string(e1) + ":  eta=" + string(factor) +"\n");
        else
            obj.physics.Update(dx);
        end
        
        % re-assemble system for new state
        obj.physics.Assemble();
        obj.physics.Constrain();

		% convergence criterion check
        if (En_err0 < 0)
            En_err0 = sum(abs(obj.physics.fint.*dx));
            En_err = En_err0;

            if (En_err0==0)
                En_err0 = 1e-12;
            end
        else
            En_err = sum(abs(obj.physics.fint.*dx));
        end
        En_err_n = En_err/En_err0;
		obj.convergence_log(stepnum,it+1) = En_err_n;
                    
        fprintf("    Residual:" + string(En_err_n) + "   ("+string(En_err)  +") \n");
        
		%determine whether to continue or exit solver loop
        it=it+1;
        if (it>curr_max_it || En_err_n<obj.Conv || En_err<obj.tiny)
            obj.physics.Commit("Pathdep");
            irr = obj.physics.Irreversibles();
            if (irr == false)
                stop = true;
            else
                obj.physics.Assemble();
                obj.physics.Constrain();
                recalc_pre = true;
                En_err0 = -1;
                curr_max_it = it + obj.maxIt;
            end
        end
    end
    
	%commit time step state and history-dependent variables
    obj.physics.Commit("Timedep");
end

