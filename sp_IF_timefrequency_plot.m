% References:   
% (1) Breaking the Limits – Redefining the Instantaneous Frequency, http://arxiv.org/abs/1605.00975
% (2) The Fourier Decomposition Method for nonlinear and nonstationary time
% series analysis, http://arxiv.org/abs/1503.06675

function [freq, amp]=sp_IF_timefrequency_plot(x,t,Fs,fmax)     
    [m,n] = size(x); 
    if(m<n) % convert row data to coloumn data
       x=x';
    end
    
    [m1,n1] = size(t);
    if(m1<n1) % convert row time to coloumn time
       t=t';
    end
    [m,n] = size(x);
    
    frequency=zeros(m,n); % initialize variables
    amplitude=frequency;  % initialize variables
    
    for i=1:n % compute freq and amp for all bands or signals in colomn
       [amp freq]=computeAmpFreq(x(:,i),Fs,t);
       amplitude(:,i)=amp;
       frequency(:,i)=freq;       
    end
              
    sp_PlotTF_frqAmp(frequency,amplitude,t,Fs,fmax);   
    
function [amp freq]=computeAmpFreq(x,Fs,t)
    plotFig=0;
    if 1    
        mx=0; mx=mean(x); % calculate mean    
        x=x-mx; % remove mean
        z=hilbert(x);
    else
        z = sp_IDCT_FSAS(x);
    end    
    amp=abs(z);
    phi=unwrap(angle(z));   
    
    if 1 % forward diff is better      
        diffPhase=diff(phi);
        if(plotFig==1)
            freq=diffPhase*(Fs/(2*pi)); 
            plot(t(1:length(freq)),freq); 
            FntSize=18;
            xlabel('Time (s)')
            ylabel('Frequency (Hz)')
            set(gca,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'xlabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'ylabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
        end
        % if diff is -ve, make it +ve by adding pi
        index=find(diffPhase<0);
        diffPhase(index)=diffPhase(index)+pi;
        diffPhase=[diffPhase;diffPhase(end)];
        if(plotFig==1)
            freq=diffPhase*(Fs/(2*pi)); 
            plot(t(1:length(freq)),freq); 
            FntSize=18;
            xlabel('Time (s)')
            ylabel('Frequency (Hz)')
            set(gca,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'xlabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'ylabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
        end        
    else % CentralDiff        
        tmp=phi;
        diffPhase=((tmp(3:end)-tmp(1:end-2))/2);
        %diffPhase=[tmp(1) tmp tmp(end)];
        if(plotFig==1)
            freq=diffPhase*(Fs/(2*pi)); 
            plot(t(1:length(freq)),freq); 
        end
        % if diff is -ve, make it +ve by adding pi
        index=find(diffPhase<0); 
        diffPhase(index)=diffPhase(index)+pi;        
        diffPhase=[diffPhase(1);diffPhase;diffPhase(end)];
        if(plotFig==1)
            freq=diffPhase*(Fs/(2*pi)); 
            plot(t(1:length(freq)),freq); 
        end
    end
    freq=diffPhase*(Fs/(2*pi));    