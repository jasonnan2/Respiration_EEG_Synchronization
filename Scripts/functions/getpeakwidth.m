function [wxPk] = getpeakwidth(y,x,iPk,wRef)
  iLBh=[1; iPk(1:end-1)];
  iRBh=[iPk(2:end); length(x)];
  % get the width boundaries of each peak
  wxPk = getUserBounds(y, x, iPk, iLBh, iRBh, wRef);
end
%--------------------------------------------------------------------------
function bounds = getUserBounds(y, x, iPk, iLB, iRB, refHeight)
  n = length(iPk);
  bounds = zeros(n,2,'like',x);
  % interpolate both the left and right bounds clamping at borders
  for i=1:n
    base(i)=refHeight;
    % compute the index of the left-intercept at half max
    iLeft = findLeftIntercept(y, iPk(i), iLB(i), refHeight);
    if iLeft < iLB(i)
      bounds(i,1) = x(iLB(i));
    else
      bounds(i,1) = linterp(x(iLeft),x(iLeft+1),y(iLeft),y(iLeft+1),y(iPk(i)),base(i));
    end
    % compute the index of the right-intercept
    iRight = findRightIntercept(y, iPk(i), iRB(i), refHeight);
    if iRight > iRB(i)
      bounds(i,2) = x(iRB(i));
    else
      bounds(i,2) = linterp(x(iRight), x(iRight-1), y(iRight), y(iRight-1), y(iPk(i)),base(i));
    end
  end
end

%--------------------------------------------------------------------------
function idx = findLeftIntercept(y, idx, borderIdx, refHeight)
  % decrement index until you pass under the reference height or pass the
  % index of the left border, whichever comes first
  while idx>=borderIdx && y(idx) > refHeight
    idx = idx - 1;
  end
end

%--------------------------------------------------------------------------
function idx = findRightIntercept(y, idx, borderIdx, refHeight)
  % increment index until you pass under the reference height or pass the
  % index of the right border, whichever comes first
  while idx<=borderIdx && y(idx) > refHeight
    idx = idx + 1;
  end
end

%--------------------------------------------------------------------------
function xc = linterp(xa,xb,ya,yb,yc,bc)
  % interpolate between points (xa,ya) and (xb,yb) to find (xc, 0.5*(yc-yc)).
  xc = xa + (xb-xa) .* (bc-ya) ./ (yb-ya);

  % invoke L'Hospital's rule when -Inf is encountered. 
  if isnumeric(xc) && isnan(xc) || coder.target('MATLAB') && isdatetime(xc) && isnat(xc)
    % yc and yb are guaranteed to be finite. 
    if isinf(bc)
      % both ya and bc are -Inf.
      if isnumeric(xa)
        xc(1) = 0.5*(xa+xb);
      else
        xc(1) = xa+0.5*(xb-xa);
      end
    else
      % only ya is -Inf.
      xc(1) = xb;
    end
  end
end