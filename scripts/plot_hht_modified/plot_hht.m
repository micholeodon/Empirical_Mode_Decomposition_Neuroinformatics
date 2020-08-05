function plot_hht(imf,Ts)
% Plot the HHT.
% plot_hht(x,Ts)
% 
% :: Syntax
%    The array x is the input signal and Ts is the sampling period.
%    Example on use: [x,Fs] = wavread('Hum.wav');
%                    plot_hht(x(1:6000),1/Fs);
% Func : emd
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


% Get HHT.
for k = 1:size(imf,1)
   b(k) = sum(imf(k,:).*imf(k,:));
   th   = angle(hilbert(imf(k,:)));
   d{k} = diff(th)/Ts/(2*pi);
end
[u,v] = sort(-b);
b     = 1-b/max(b);

% Set time-frequency plots.
N = size(imf,2);
c = linspace(0,(N-2)*Ts,N-1);
for k = v(1:2)
   figure, plot(c,d{k},'k.','Color',b([k k k]),'MarkerSize',3);
   set(gca,'FontSize',8,'XLim',[0 c(end)],'YLim',[0 1/2/Ts]); xlabel('Time'), ylabel('Frequency');
end
