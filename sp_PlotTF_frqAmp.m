function [nt,tscale,fscale]=sp_PlotTF_frqAmp(freq,amp,t,Fs,fmax)
    
    %t0=0;t1=t(end); 
    t0=t(1);t1=t(end);
    multfactor=6;
    if(length(freq(:,1))>=100*multfactor) 
        fres=100*multfactor; tres=100*multfactor;
    else
        fres=length(freq(:,1));    tres=fres;
    end
    tres=length(t); % ridge extraction
    %fres=length(freq(:,1));    tres=fres;
    
    %fres=268;    tres=fres; % for earthquake data
    
    %fw0=0; fw1=max(max(freq));
    if 0
        fw0=min(min(freq)); fw1=max(max(freq));
        if(fw0<0)          fw0=0;    end
    else % Graviattion wave data
        fw0=min(min(freq)); fw1=fmax;
        if(fw0<0)          fw0=0;    end    
    end
    %fw1=Fs/2; % max frequency in plot
    tw0=t0;     tw1=t1;
    %  4.call nspplote.m to plot the figure of time-frequency spectrum 
    %----- Get the values to plot the spectrum
    lscale=0;
    [nt,tscale,fscale] = nspplote(freq,amp,t0,t1,fres,tres,fw0,fw1,tw0,tw1,lscale); 
    %for clear appearance,clear the low-freqency if you want 
    setclear=0; 
    if setclear==1
        nt(1,:)=0;
    end
    
    q=fspecial('gaussian',7,0.6);  
    nsu=filter2(q,nt);
    nsu=filter2(q,nsu);
    figure;
    %subplot(5,1,1);
    imagesc(tscale,fscale,nsu.^.5); 
    axis xy;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    %title('FMD based Time Frequncy Plot');
    FntSize=18;
    title('Time-Frequency-Energy plot','FontSize',FntSize,'FontName','Times'); 
    %% Use end code of RCADA
    
    
    set(gca,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'xlabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'ylabel');
    set(h1,'FontSize',FntSize,'FontName','Times') 
    
    if 0% marginal spectrum plot
        %marginal spectrum
        ms=sum(nt,2);
        hsp_fre1=fw1; %give frequency axis max value:
        %comparisoncoefficients
        hco=Fs/2/hsp_fre1;
        figure;
        plot(fscale(1:end),hco*ms(1:end));
        %loglog(fscale(1:end),hco*ms(1:end));
        xlabel('Frequency (Hz)')
        ylabel('Spectral Density')
        title('Marginal spectrum','FontSize',16,'FontName','Times');    
        FntSize=16;
        set(gca,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'xlabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'ylabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
    
    end
    if 0% E(t) plot, Energy fluctuation with time
        %marginal spectrum
        ms=sum(nt,1);
        hsp_fre1=fw1; %give frequency axis max value:
        %comparisoncoefficients
        hco=Fs/2/hsp_fre1;
        figure;
        plot(tscale,hco*ms(1:end));
        %loglog(fscale(1:end),hco*ms(1:end));
        xlabel('Time (s)')
        ylabel('Energy Density')
        title('Instantaneous energy E(t)','FontSize',FntSize,'FontName','Times');    
        %FntSize=16;
        set(gca,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'xlabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'ylabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
    
    end