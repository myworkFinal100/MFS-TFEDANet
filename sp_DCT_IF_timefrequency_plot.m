% References:   
% (1) Breaking the Limits ? Redefining the Instantaneous Frequency, http://arxiv.org/abs/1605.00975
% (2) The Fourier Decomposition Method for nonlinear and nonstationary time
% series analysis, http://arxiv.org/abs/1503.06675

function [frequency, amplitude]=sp_DCT_IF_timefrequency_plot(x,t,Fs,fmax)      
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
    mx=0; mx=mean(x); % calculate mean    
    x=x-mx; % remove mean to obtain meaningfull IF, e.g., see delta[n-n)] example
    if 0 % Hilbert based Analytic signal cmputation 
        z=hilbert(x);
    else % FSAS
        if 0 % direct computation, which is computationally complex
            %%%%%%%%%%%% DCT based Analytic signal cmputation
            y=dct(x);
            z=0;
            N=length(x);
            n=1:N;
            w=ones(1,N)*(2/N)^(0.5); w(1)= (1/N)^(0.5);       
            for k=1:N
                z= z+w(k)*y(k)*exp(1i*pi*(2*n-1)*(k-1)/(2*N) );
            end
            z=z.';
            %%%%%%%%%%%%
        else % FFT based FSAS computation
            z = sp_IDCT_FSAS(x);
        end
    end
    amp=abs(z);
    if 0
        phi=unwrap(angle(z)); 
    else        
        phi=sp_PhaseUnwrap(angle(z)); 
    end
    
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
    freq=medianfilter(freq,5);    
    %freq=smooth(freq);
    