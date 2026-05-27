classdef CoulombFriction < handle
  
  
  properties
    frictionAngle
    cohesion
    elasticSlip             % state variable for computing traction with penalty
    stiffN
    stiffT
  end
  
  methods
    function obj = CoulombFriction(input)

      default = struct('cohesion',0.0,'frictionAngle',30.0,...
                       'normalStiff',kN,'tangentialStiff',kT);

      params = readInput(input);

      obj.cohesion = params.cohesion;
      obj.frictionAngle = deg2rad(params.frictionAngle);
      obj.elasticSlip = 3*size(obj.cohesion,1);

    end

    
    function [t,dtdg] = updateTraction(obj,f,gCurr,gOld,state)
      % update traction and compute derivative w.r.t gap

      % tOld: last converged traction
      % jumpNew/jumpOld: current and last converged fracture jump

      t = zeros(3,1);
      dtdg = zeros(3);

      if state == ContactMode.open
        % open dof stay open
        return
      end

      t(1) = obj.stiffN * gCurr(1);

      slip = gN(2:3) - gOld(2:3);
      elasticSlip = obj.elasticSlip([2*f-1;2*f]);

      % trial traction
      tT = [obj.penalty_t * elasticSlip + slip];


      tractionNew(1) = tTrial(1);
      dtdg(1,1) = obj.stiffN;

      tauLim = obj.cohesion(f) - tan(obj.frictionAngle(f))*t(1);


      if state == ContactMode.stick

        t(2:3) = tT;
        dtdg(2,2) = obj.stiffT;
        dtdg(3,3) = obj.stiffT;


      elseif state == ContactMode.slip || state == ContactMode.newSlip

        % slipNorm is always 0 at first iteration
        slipNorm = norm(slip);
        isSlidingReliable = slipNorm > obj.activeSet.tol.sliding;

        dtdg(1,1) = obj.penalty_n;

        if isSlidingReliable % we can trust the available tangential traction direction

          slipDir = tTrial_t / tTrial_t_norm; 

          % consistent tangent operator
          dtdg(2:3,2:3) = obj.penalty_t * tauLim * (tTrial_t_norm^2*eye(2) - tTrial_t * tTrial_t')/tTrial_t_norm^3;

          % slipDir = slip/slipNorm;
          % dtdg(2:3,2:3) = tauLim * (slipNorm^2*eye(2) - slip * slip')/slipNorm^3;
        
        
        
        else

          slipDir = tCurr(2:3)/norm(tCurr(2:3));
          
        end

        tractionNew(2:3) = tauLim * slipDir;
        dtdg(2:3,1) = -obj.penalty_n * tan(obj.phi(fracId)) * slipDir;

      end

    end
  end
end

