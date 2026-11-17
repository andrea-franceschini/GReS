function fillStats(obj,Tend,time,symValue,T_setup,T_setup_eff,setupTSAM,tolReq,tolScaled,nonlinIter,tolReach)
% FILLSTATS Records timing, iteration count, and residual statistics on linearSolver obj

   obj.aTimeSolve = obj.aTimeSolve + Tend;
   obj.nSolve = obj.nSolve + 1;
   solves = obj.nSolve;
   obj.aIter = obj.aIter + obj.params.iter;
   obj.maxIter = max(obj.maxIter,obj.params.iter);
   obj.precCompLin(solves) = T_setup;
   obj.precEffLin(solves) = T_setup_eff;
   obj.SAM.CompLin(solves) = setupTSAM;
   obj.iterLin(solves) = obj.params.iter;
   obj.timeLin(solves) = time; 
   obj.solveTLin(solves) = Tend;
   obj.symFlagLin(solves) = symValue;
   obj.tolReqLin(solves) = tolReq;
   obj.scaledTolLin(solves) = tolScaled;
   obj.newtonLin(solves) = nonlinIter;
   obj.tolReachLin(solves) = tolReach;
end
