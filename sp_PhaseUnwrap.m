% slightly different than matlab implementation and One-Dimensional Phase Unwrapping Problem By Dr. Munther Gdeisat and Dr. Francis Lilley
function xu=sp_PhaseUnwrap(xw)    
    xu = xw;           
        tol1=10^(-10); % this gives better performance  
        %tol2=10^(-12);
        for i=2:length(xw)
            difference = xw(i)-xw(i-1);
            if( (difference > (pi+tol1)) )
                xu(i:end) = xu(i:end) - 2*pi;
            elseif( (difference <= -(pi-tol1)) )
                xu(i:end) = xu(i:end) + 2*pi;
            end
        end
    