function imf = emd(x)
% Empiricial Mode Decomposition (Hilbert-Huang Transform)
% imf = emd(x)
% Func : findpeaks
%
% Based on: Alan Tan, 2006
% Modified  02.04.2020 by: Michal Konrad Komorowski, michu.kom@gmail.com
%
%
% Copyright (c) 2016, Alan Tan
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

imf = [];
while ~ismonotonic(x)
    x1 = x;
    sd = Inf;
    nSifts = 0;
    while (sd > 0.1) | ~isimf(x1)
        nSifts = nSifts + 1;
        s1 = getspline(x1);
        s2 = -getspline(-x1);
        x2 = x1-(s1+s2)/2;
        
        sd = sum((x1-x2).^2)/sum(x1.^2);
        x1 = x2;
    end
    
%     imf{end+1} = x1;
    imf = [imf;x1];
    x          = x-x1;
end
% imf{end+1} = x;
imf = [imf;x];


% FUNCTIONS

    function isMonotonic = ismonotonic(x)
        
        % will equal zero if only one type of extrema is found
        test = length(findpeaks(x))*length(findpeaks(-x)); 
        
        if test > 0
            isMonotonic = 0;
        else
            isMonotonic = 1;
        end
    end

    function isIMF = isimf(x)
        
        N  = length(x);
        % 1. The number of extrema and the number of zero crossings must either equal
        % or differ at most by one
        nZeroCrss = sum(x(1:N-1).*x(2:N) < 0); % number of zero-crossings
        nExtr = length(findpeaks(x))+length(findpeaks(-x)); % number of extrema
        
        cond1 = (abs(nZeroCrss-nExtr) <= 1);
        
        % 2. At any point, the mean value of the envelope defined by the 
        % local maxima and the envelope defined by the local minima
        % is zero
        %
        % In practice this condition is too strong, so it is subsituted by SD criterion presented in Huang et al. 1996 eq. 5.5
        % which is present in the main function
        [emin, emax] = getenvelopes(x);
        tol = 1e-7;
        cond2 = all( ((emax+emin)/2) < tol);
        
        if cond1 % && cond2
            isIMF = 1;
        else
            isIMF = 0;
        end
    end

    function s = getspline(x)
        
        N = length(x);
        p = findpeaks(x);
        s = spline([0 p N+1],[0 x(p) 0],1:N);
    end

    function [emin, emax] = getenvelopes(x)
        
        N = length(x);
        p = findpeaks(x);
        emax = spline([0 p N+1],[0 x(p) 0],1:N);
        p = findpeaks(-x);
        emin = spline([0 p N+1],[0 x(p) 0],1:N);
    end

end