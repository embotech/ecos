function [cnew,Gnew,hnew,dimsnew,Anew] = prob2fixedcones(c,G,h,dims,A,CONESIZE)
% Recursive function that maps an SOCP of arbitrary cone sizes into an SOCP
% with a fixed conesize, possibly increasing the number of cones.

% do nothing if the problem does not have second-order cones
if( isempty(dims.q) )
    return;
end

Gnew = G(1:dims.l,:);
hnew = h(1:dims.l,:);
n=size(G,2);
dimsnew.l = dims.l;

newvar = 0;

% Go about this cone by cone. Each CONESIZE variables are put into a new
% cone.
k = 1;
for i = 1:length(dims.q)
    
    coneidx = dims.l + sum(dims.q(1:i-1)) + (1:dims.q(i));
    Gi = G(coneidx,:);
    hi = h(coneidx,:);
    thisconesize = dims.q(i);
    
    % CASE 1: current cone size is smaller than CONESIZE - introduce new
    % fake variables.
    if( thisconesize < CONESIZE )
        padding = CONESIZE - thisconesize;
        Gnew = [Gnew;
            Gi, zeros(thisconesize,newvar);
            zeros(padding,n+newvar)];
        hnew = [hnew; hi; zeros(padding,1)];
        dimsnew.q(k) = CONESIZE;
        k = k+1;
    end
    
    % CASE 2: current cone size equals CONESIZE - just copy new cone in
    if( thisconesize == CONESIZE )
        Gnew = [Gnew;
            Gi, zeros(thisconesize,newvar)];
        hnew = [hnew; hi];
        dimsnew.q(k) = CONESIZE;
        k = k+1;
    end
    
    % CASE 3: current cone size is larger than CONESEIZE - split down into
    % new cones of CONESIZE or less
    if( thisconesize > CONESIZE )
        
        % generate first new cones, and put the remainder at the bottom
        newcones = floor((thisconesize-1)/(CONESIZE-1));
        for j = 1:newcones
            newvar = newvar+1;
            newconeidx = 1 + (j-1)*(CONESIZE-1) + (1:(CONESIZE-1));
            Gnew = [Gnew, zeros(size(Gnew,1),1);
                zeros(1,size(Gnew,2)) , -1;
                Gi(newconeidx,:), zeros(CONESIZE-1,newvar) ]; %#ok<*AGROW>
            hnew = [hnew; 0; hi(newconeidx)];
            dimsnew.q(k) = CONESIZE;
            k = k+1;
        end
        
        % a new cone gets anything that didn't get cast into new ones
        Gnew = [Gnew;
            Gi(1,:), zeros(1,newvar);
            zeros(newcones,n+newvar-newcones), -eye(newcones)];
        hnew = [hnew; hi(1); zeros(newcones,1)];
        dimsnew.q(k) = 1 + newcones;
        
        remainder = newconeidx(end)+1:thisconesize;
        if( ~isempty(remainder) )
            Gnew = [Gnew; Gi(remainder,:), zeros(length(remainder),newvar)];
            hnew = [hnew; hi(remainder)];
            dimsnew.q(k) = dimsnew.q(k) + length(remainder);
        end
        k = k+1;
    end
end

Anew = [A, zeros(size(A,1),newvar)];
cnew = [c; zeros(newvar,1)];


assert(length(hnew)==dimsnew.l+sum(dimsnew.q),'h and dims mismatch');
assert(size(Gnew,1)==dimsnew.l+sum(dimsnew.q),'G,1 and dims mismatch');
assert(size(Gnew,2)==n+newvar,'G,2 and nvar mismatch');
assert(size(Anew,2)==n+newvar,'A,2 and nvar mismatch');
assert(size(cnew,1)==n+newvar,'c,1 and nvar mismatch');

% recursive call
if( any(dims.q ~= CONESIZE) )
    [cnew,Gnew,hnew,dimsnew,Anew] = prob2fixedcones(cnew,Gnew,hnew,dimsnew,Anew,CONESIZE);
end

