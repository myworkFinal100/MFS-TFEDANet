
function tt=sp_GravitationalWaveAnalysis1()
    close all; 
    clear all;
    %format shortEng;
    clc
    
    if 0 % for oscilloscope
        Fs=48000;
        t=0:1/Fs:10-1/Fs;
        f0=20; t1=t(end); f1=10000;
        y = chirp(t,f0,t1,f1);
        plot(t,y);
        axis([t(1) t(end) -1.5 1.5]);
        variables(:,1)=t;
        variables(:,2)=y;
        save('Osc.mat','variables')
        
        %yy=spectrogram(y,256,250,256,Fs,'yaxis');        
        [~,~,~,pxx,fc,tc] = spectrogram(y,kaiser(128,18),120,128,Fs, ...
        'MinThreshold',-30);
        plot(tc(pxx>0),fc(pxx>0),'.');        
        pp=1;
    end
    
    G=(6.67191)*10^(-11); 
    Mo=(1.98855)*10^(30); % solar mass (M?)
    %M=62*Mo;
    m1=36.2*Mo;
    m2=29.1*Mo;
    M=m1+m2;
    c=299792458;    
    
    if 1
        x1=load('J:\JIIT2BU\PapersToBeSubmitted\GavitationalWaveTF_withAmit\GravitationalWave\data\fig1-observed-H_FLD.txt');
        x1nr=load('J:\JIIT2BU\PapersToBeSubmitted\GavitationalWaveTF_withAmit\GravitationalWave\data\fig1-waveform-H_NR_FLD.txt');
    else
        x1=load('J:\JIIT2BU\PapersToBeSubmitted\GavitationalWaveTF_withAmit\GravitationalWave\data\fig1-observed-L_FLD.txt');
        x1nr=load('J:\JIIT2BU\PapersToBeSubmitted\GavitationalWaveTF_withAmit\GravitationalWave\data\fig1-waveform-L_NR_FLD.txt');
    end
    %x1=load('J:\JIIT2BU\PapersToBeSubmitted\GavitationalWaveTF_withAmit\GravitationalWave\data\fig2-unfiltered-waveform-H_NR_FLD.txt');

    % second event data
    %x1=load('J:\JIIT2BU\PapersToBeSubmitted\GavitationalWaveTF_withAmit\GravitationalWave\data1\GW150914_tutorial\GW150914_4_NR_waveform.txt');

    dsize=min(size(x1(:,1)),size(x1nr(:,1)));
    x1=x1(1:dsize,:); x1nr=x1nr(1:dsize,:);

    t=x1(:,1);
    x=x1(:,2);
    
    %Horg=stat41(x');
    
    [maxdataValue MaxIndex]=max(abs(x));  
    maxt=t(MaxIndex);
    %x=x-mean(x);% remove mean   
    Fs=1/(t(2)-t(1)); 
    tt=1;
    L=length(x);
    %NFFT = 2^nextpow2(L); % Next power of 2 from length of y 
    NFFT=1024*32*8;
    Y = fft((x),NFFT)/L;           
    f = Fs/2*linspace(0,1,NFFT/2+1);  
    ampf=2*abs(Y(1:NFFT/2+1));
    %plot(f,2*abs(Y(1:NFFT/2+1)));
    %plot(amp(1:130),f(1:130));
    range=1000*8;         
    fmax=400;       

    if 1 % wsst
        [sst,f1] = wsst(x,Fs);
        %[fridge,iridge] = wsstridge(sst);
        [fridge,iridge] = wsstridge(sst,10,f1,'NumRidges',2);
        %pcolor(t,f,abs(sst));        shading interp
        contour(t,f1,abs(sst));
        title('Synchrosqueezed Transform')
        hold on
        %plot(t,fridge)
        hold on;
        plot(t,fridge,'k','linewidth',2);
        title('Synchrosqueezed Transform with Overlaid Ridge')
        tt=1;  
    end
    
    
    if 1          
         figure
         subplot(2,1,1);
         plot(t,x); 
         FntSize=18; 
         axis([t(1) t(end) min(x) max(x)])
         grid on;
         ylabel('Strain (10^{-21})')
         xlabel('Time (s)')
         set(gca,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'xlabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'ylabel');
         set(h1,'FontSize',FntSize,'FontName','Times')

         subplot(2,1,2);
         %figure
         %plot(f(1:range),ampf(1:range));

         %plot(f(1:range),20*log(ampf(1:range)));
         plot(f(1:range),ampf(1:range));
         %axis([f(1) f(end) -1.5 1.5])
         grid on;
         ylabel('|H1(f)|(10^{-21})')
         xlabel('Frequency (Hz)')
         set(gca,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'xlabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'ylabel');
         set(h1,'FontSize',FntSize,'FontName','Times')   
         %subplot(3,1,3);

    end
    
    if 1 % TFE by orthogonal (dft) or LINOEP (fir) decomposition                    
       CutFrq1=[350 300 200 100 60 25]; 
       ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(x,t,Fs,CutFrq1,'dct');
       figure;
       for i=1:length(CutFrq1)+1 
           Data(:,i)=ImfByFilter{i}; 
           subplot(length(CutFrq1)+1,1,i);
           plot(t,Data(:,i));
           tmp=int2str(7-i);
           tmp1=strcat('FIBF',tmp);
           ylabel(tmp1)
           set(gca,'FontSize',FntSize,'FontName','Times')
           h1 = get(gca, 'xlabel');
           set(h1,'FontSize',FntSize,'FontName','Times')
           h1 = get(gca, 'ylabel');
           set(h1,'FontSize',FntSize,'FontName','Times')
           axis([t(1) t(end) min(Data(:,i)) max(Data(:,i))])
           if(i<length(CutFrq1)+1)
               set(gca,'xticklabel',{[]}) 
           end
       end
    end
    Data1=Data(:,2:end-1); 
    Len=length(Data1(1,:));
    figure;
    for i=1:Len
        subplot(Len+1,1,i);
        data=Data1(:,i);
        plot(data);
        hold on
        [val ind]=max(abs(data)); 
        mu=t(ind); 
         
         %sigma=var(data); % OK
         sigma=0.005*i; % OK
         
        win=exp(-(t-mu).^2/(2*sigma^2));
        if(i==1)
            win = sp_gaussianwin(length(t),.25,0.01); % (n,r,sigma)
            win(1:end/2)=1;
            win=1-win;
         end
        plot(win);
        if(i==Len)
            [val, ind]=max((win)); 
            win(1:ind)=val;
            plot(win);
        end
        data=data.*win;
        hold off
        plot(t,data);
        set(gca,'FontSize',FntSize,'FontName','Times')
           h1 = get(gca, 'xlabel');
           set(h1,'FontSize',FntSize,'FontName','Times')
           h1 = get(gca, 'ylabel');
           set(h1,'FontSize',FntSize,'FontName','Times')
           axis([t(1) t(end) min(data) max(data)])
           if(i<=Len)
               set(gca,'xticklabel',{[]}) 
           end
        Data1(:,i)=data;
        pp=1;
    end
    subplot(Len+1,1,i+1);    
    gwdata=sum(Data1,2);
    plot(t,gwdata);
    axis([t(1) t(end) min(gwdata) max(gwdata)]);
    set(gca,'FontSize',FntSize,'FontName','Times')
           h1 = get(gca, 'xlabel');
           set(h1,'FontSize',FntSize,'FontName','Times')
           h1 = get(gca, 'ylabel');
           set(h1,'FontSize',FntSize,'FontName','Times')
    
    figure; 
    fsst(gwdata,Fs,'yaxis');
    wsst(gwdata,Fs);
    subplot(2,1,1)
    plot(t,gwdata,'LineWidth',2);
    hold on
    plot(x1(:,1),x1(:,2),'r.','LineWidth',2);
    axis([t(1) t(end) -1.3 1.3])
    ylabel('Strain (10^{-21})')
    %xlabel('Time (s)')
    set(gca,'XTick',[]);
    legend('Proposed reconstruction','GW Strain time-series');
    FntSize=20;
    set(gca,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'xlabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'ylabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    
    subplot(2,1,2) 
    residue=x1(:,2)-gwdata;
    plot(t,residue,'LineWidth',2);
    %hold on
    %plot(x1nr(:,1),x1nr(:,2),'r','LineWidth',2);
    axis([t(1) t(end) -1.3 1.3])
    ylabel('Strain (10^{-21})')
    xlabel('Time (s)')
    legend('Residue');
    FntSize=20;
    set(gca,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'xlabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'ylabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    
    figure;       
    subplot(2,1,1)
    plot(t,gwdata,'LineWidth',2);
    hold on
    plot(x1nr(:,1),x1nr(:,2),'r.','LineWidth',2);
    axis([t(1) t(end) -1.3 1.3])
    ylabel('Strain (10^{-21})')
    %xlabel('Time (s)')
    set(gca,'XTick',[]);
    legend('Proposed reconstruction','Numerical relativity (NR)');
    FntSize=20;
    set(gca,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'xlabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'ylabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    
    subplot(2,1,2)
    plot(t,gwdata-x1nr(:,2),'LineWidth',2);
    axis([t(1) t(end) -1.3 1.3])
    ylabel('Strain (10^{-21})')
    xlabel('Time (s)')
    legend('Difference between estimated and NR');
    FntSize=20;
    set(gca,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'xlabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'ylabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    
    %plot(t,x,'k');    
    %plot(t,abs(hilbert(gwdata)))
    
    if 1 % wsst
        figure
        [sst,fr] = wsst(x,Fs);
        %[fridge,iridge] = wsstridge(sst);
        [fridge,iridge] = wsstridge(sst,10,fr,'NumRidges',3);
        %pcolor(t,f,abs(sst));        shading interp
        contour(t,fr,abs(sst));
        title('Synchrosqueezed Transform') 
        hold on
        %plot(t,fridge)
        hold on;
        plot(t,fridge,'k','linewidth',2);
        title('Synchrosqueezed Transform with Overlaid Ridge')
        xrec1 = iwsst(sst,iridge(:,1));
        xrec2 = iwsst(sst,iridge(:,2));
        figure
        plot(xrec1); hold on;
        plot(xrec2);
        xrec3=[xrec1(1:1577);xrec2(1578:2792);xrec1(2793:end)];
        
        figure
        [sst,fr] = wsst(xrec1,Fs);
        %[fridge,iridge] = wsstridge(sst);
        [fridge,iridge2] = wsstridge(sst,10,fr,'NumRidges',1);
        %pcolor(t,f,abs(sst));        shading interp
        contour(t,fr,abs(sst));
        title('Synchrosqueezed Transform')
        hold on
        %plot(t,fridge)
        hold on;
        plot(t,fridge,'k','linewidth',1);
        title('Synchrosqueezed Transform with Overlaid Ridge')
        
        xrec12 = iwsst(sst,iridge2);
        figure
        plot(xrec12);
         
        tt=1;   
    end

    
    fmax=500;
    [freq1 amp1]=sp_DCT_IF_timefrequency_plot(x,t,Fs,fmax);
    [freq1 amp1]=sp_DCT_IF_timefrequency_plot(x1nr(:,2),t,Fs,fmax);
    [freq1 amp1]=sp_DCT_IF_timefrequency_plot(gwdata,t,Fs,fmax);
     %snr=70; y = awgn(gwdata,snr,'measured');
     envelop=abs(sp_IDCT_FSAS(gwdata));
    gwdata_unitAmp=gwdata./envelop; 
    [freq1 amp1]=sp_DCT_IF_timefrequency_plot(gwdata_unitAmp,t,Fs,fmax);
    wsst(gwdata_unitAmp,Fs);
    tt=1;    
    if 0
        nr=x1nr(:,2); 
        Hnr=stat41(nr');
        Hdct=stat41(gwdata');
        
        filename = 'org_nr_proposed.xlsx';
        A(:,1)=t; A(:,2)=x; A(:,3)=x1nr(:,2); A(:,4)=gwdata;
        xlswrite(filename,A)
    end
    
    
    if 1
        
        [freq amp]=sp_IF_timefrequency_plot(x1nr(:,2),t,Fs,fmax);
        [freq amp]=sp_IF_timefrequency_plot(gwdata,t,Fs,fmax); 
        snr=60;
        gwdata = awgn(gwdata,snr,'measured');
        [freq amp]=sp_IF_timefrequency_plot(gwdata,t,Fs,fmax);     
        figure; 
        plot(t, freq);   
        
        indx=find(freq>400); % value > 300 set them zero
        freq(indx)=0;        
        plot(t, freq);
        indx=find(freq(1:end/10)>40); % value > 300 set them zero 
        freq(indx)=0;
        plot(t, freq);
        
        
        start=30; stop=300;  
        freq_ss=zeros(length(t),1);
        freq_ss(1)=30;%freq_ss(end)=300;
        freq(1)=30; %freq(end)=300;        
        
        tmp=30;        
        for i=2:length(freq)-1             
            if( ( freq(i) > tmp ))
               freq_ss(i)=freq(i);
               tmp=freq(i);
            else
                freq_ss(i)=tmp;
            end
        end
        
        figure;
        
        plot(t,freq_ss);
        axis([t(1) t(end) min(freq_ss) max(freq_ss)]);
        
        indx=find(freq_ss);
        val=freq_ss(indx);
        envlp= spline(indx,val,1:length(t)); %upper spline bound of this sift
        hold on
        plot(t,envlp); 
        
       pp=1; 
    end
    
    if 0
        figure
        y=gwdata;
        spectrogram(y,kaiser(128,18),120,128,Fs,'reassigned','yaxis');
        
        [~,~,~,pxx,fc,tc] = spectrogram(y,kaiser(128,18),120,128,Fs, ...
        'MinThreshold',-30);
        plot(tc(pxx>0),fc(pxx>0),'--');
        pp=1; 
    end
    
    if 1 
    figure;   
    plot(t,gwdata);
    hold on; 
    dd=1:1:length(t);
    [spmax, spmin, flag]=extrema(gwdata);  %call function extrema           
    upper= spline(spmax(:,1),spmax(:,2),dd); %upper spline bound of this sift 
    plot(t,upper);
    lower= spline(spmin(:,1),spmin(:,2),dd); %lower spline bound of this sift
    plot(t,lower);
    
    plot(t,(lower+upper)/2);
    
    
    %gwdata1=gwdata./(-lower'+upper')/2;
    figure; 
    gwdata1=gwdata-(lower'+upper')/2;
    plot(t,gwdata1);
    hold on
    plot(t,gwdata);
    
    figure;  
    plot(t,gwdata1);
    [spmax, spmin, flag]=extrema(gwdata1);  %call function extrema            
    upper= spline(spmax(:,1),spmax(:,2),dd); %upper spline bound of this sift 
    hold on
    plot(t,upper);
    lower= spline(spmin(:,1),spmin(:,2),dd); %lower spline bound of this sift
    plot(t,lower);
    
    [maxmin, flag]=sp_extrema(gwdata);
    
    figure
    plot(maxmin(:,1),maxmin(:,2));
    
    hold on
    
    CutFrq1=[260 30]; 
       ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(x,t,Fs,CutFrq1,'dft');
       for i=1:length(CutFrq1)+1 
           Data(:,i)=ImfByFilter{i}; 
           subplot(length(CutFrq1)+1,1,i);
           plot(Data(:,i));
       end
       
       filtered_x=Data(:,2);
       
    figure
    dataloc=maxmin(:,1);
    dataval=filtered_x(dataloc);
    plot(dataloc,dataval);
    splineint= spline(dataloc,dataval,dd);
    hold on
    plot(splineint);
    
    figure
    plot(gwdata-x);
    hold on
    plot(gwdata); 
    
    pp=1;   
    end
        
    
    [freq amp]=sp_IF_timefrequency_plot(gwdata,t,Fs,fmax);    
    figure;
    plot(t, freq);
    
    sp_TimeFreqPlotInputIMFsFIBFs(x,t,fmax) ; % pass colomn vectors, emd method of IF calculation   

    if 1
        %----- Get frequency and amplitude
        dt=t(2)-t(1);ifmethod=[];normmethod=[];nfilter=[];
        [freq1,amp1] = sp_fa(gwdata,dt,ifmethod,normmethod,nfilter);
        plot(t,freq1); 
    end
    if 1 
        [freq amp]=sp_IF_timefrequency_plot(x,t,Fs,fmax);  % add data to combat end effects        
        %[minFreqValue MinIndex]=min(freq(1:MaxIndex));
        %freq(1:MinIndex)=minFreqValue; 
    end

    if 1
        plot(t,freq);
        hold on
        plot(t,freq1,'r');
        tt=1;       
    end 

    
        tt=1  
        
      
        
    function tp = cwt_sp(x,t,Fs)
        %if 1 % wavelet analysis
            X=x;
                wname = 'morl';   
                % Define scales.
                amax = 13; scales = (1.625).^[1:0.5:amax]; % choose min scale (max frequnecy) such that its (Fs/2) and max scale (min frequency) to desired min frequency
                %scales = 1:4:128; 
                coefs = cwt(X,scales,wname,'lvlabs');  
                %freq = scal2frq(scales,wname,1/Fs);
                cfreq = centfrq(wname,16,'plot'); freq = cfreq./(scales.*(1/Fs)); % choose min scale (max frequnecy) such that its (Fs/2) and max scale (min frequency) to desired min frequency 
                figure; coefsSquared = abs(coefs).^2; imagesc(t,freq,coefsSquared); grid off; axis xy;
                figure
                surf(t,freq,abs(coefs));shading('interp'); axis tight; xlabel('Time (s)'); ylabel('Pseudo-Frequency (Hz)') 

                SCImg = wscalogram('image',(coefs),'scales',(freq),'ydata',X,'xdata',t);

                figure
                %SCImg = wscalogram('image',flipud(coefs),'scales',fliplr(freq),'ydata',X,'xdata',t);  axis tight; xlabel('Seconds'); ylabel('Pseudo-Frequency (Hz)') 
                SCImg = wscalogram('image',flipud(coefs),'scales',fliplr(freq),'ydata',X,'xdata',t);  axis tight; xlabel('Seconds'); ylabel('Pseudo-Frequency (Hz)') 
                %SCImg = wscalogram('image',flipud(coefs),'scales',fliplr(freq));  axis tight; xlabel('Seconds'); ylabel('Pseudo-Frequency (Hz)') 
            
            FntSize=18;
            set(gca,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'xlabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'ylabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            tt=1; 
        %end
        
        
            
    