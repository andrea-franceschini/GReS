
function outState = updateContactState(inState,traction,tLimit,normalGap,tols)
% tols: structure with tolerances named
% according to initializeActiveSet method
% the method is static to be reused by other contact solvERS
outState = inState;

tol = 1e-8;

% contact state update
if inState == ContactMode.open
  % from open to stick
  if normalGap < -  tols.normalGap
    outState = ContactMode.stick;
  end

elseif traction(1) > tols.normalTraction
  outState = ContactMode.open;

else % not open

  % tangential traction norm
  tau = norm(traction(2:3));
  % limiting traction

  % set to 0 a tLimit that is too small
  if tLimit < tols.minLimitTraction
    tLimit = 0;
  end

  % relax stick/sliding transition if the state is about to change
  if inState == ContactMode.stick && tau >= tLimit

    % reduce the tau if goes above limit
    tau = tau*(1-tols.tangentialViolation)-tol;

  elseif inState ~= ContactMode.stick  && tau <=tLimit

    % increase tau if falls below limit
    tau = tau*(1+tols.tangentialViolation)+tol;
  end

  % change the state after relaxation
  if tau > tLimit
    if inState == ContactMode.stick
      outState = ContactMode.newSlip;
    else
      outState = ContactMode.slip;
    end
  else
    outState = ContactMode.stick;
  end
end
end

