    function tt=sp_Speech_DCT()
    close all; 
    clear all;
    format shortEng;
    %format long;
    clc
    
    if 1
        dt=1;%365*24*60*60;
        Fs=1/dt;
        tmp1=Fs/2;
        
        if 0
                    Npart=10; % equal parts
                    inc=(Fs/2)/Npart;
                    for i=1:Npart-1
                        CutFrq1(i)=Fs/2-i*inc;
                    end
        end
        
        if 1
            for i=1:5
                        %CutFrq1(i)=(Fs/2^(i+1)); % dyadic
                        CutFrq1(i)=(Fs/2^(i+1)); % dyadic
            end
        end
        
        tmp=(1./CutFrq1)/dt;        
        tt=1; 
    end
    
    if 0
        N=2000;
        xn=zeros(1,N);
        xn(N/2)=1;
        plot(xn);
        xn=xn-mean(xn);
        plot(xn);
        Xk=dct(xn); 
        plot((Xk));
        tt=1;
    end
    
    if 1        
            %Fs=50; % for earthquake signal            
            load elcentro_NS.dat;
            x=elcentro_NS(:,2);
            %t =0:length(X)/Fs;
            t=elcentro_NS(:,1);
            Fs=1/(t(2)-t(1));
            
            fmax=Fs/2; 
            [freq1 amp1]=sp_DCT_IF_timefrequency_plot(x,t,Fs,fmax); 
        
            
            subplot(2,1,1) 
            FntSize=18;
            xlabel('Time (s)','fontsize',FntSize)
            ylabel('Acceleration (G)','fontsize',FntSize)
            axis([min(t) max(t) min(x) max(x)]);
            set(gca,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'xlabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'ylabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            
            plot(t,x);
            tt=1;
            
            L=length(x);NFFT=2^nextpow2(L);
            Y = fft(x,NFFT)/L;           
            f = Fs/2*linspace(0,1,NFFT/2+1);  
            ampf=2*abs(Y(1:(NFFT/2)+1)); 
            subplot(2,1,2);
            %SigLen=length(f();
            
            if 1 % equal energy bands
                amp2=ampf.*ampf;
                TotalEnrg=sum(amp2);
                NbOfEqEnrgBands=10;
                EnrgyInOneBand=TotalEnrg/NbOfEqEnrgBands;
                j=1;fromk=1; SigFreqLen=length(amp2);
                for i=1:SigFreqLen
                    EachBandEnergy(j)=sum(amp2(fromk:i));
                    if( (EachBandEnergy(j)>=EnrgyInOneBand) || ( abs(EachBandEnergy(j)-EnrgyInOneBand)<=EnrgyInOneBand*0.01) )
                        freqIndex(j)=i;
                        fromk=i+1;
                        j=j+1;                        
                    end
                    freqIndex(j)=SigFreqLen;
                end
                FreqHz=f(freqIndex);
                FreqHz=fliplr(FreqHz); 
                CutFrq1=0; CutFrq1=FreqHz(2:end);
            end

            [pks,locs] = findpeaks(ampf,'MinPeakDistance',40);
            pkfrq=f(locs);
            plot(f,ampf,f(locs),pks,'or')
            
            plot(f,ampf);   
            %plot(amp(1:130),f(1:130));
            xlabel('f (Hz)')
            ylabel('Speech Spectrum'); FntSize=18;
            set(gca,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'xlabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'ylabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            axis([1 max(f) min(ampf) max(ampf)]);
        
    tt=1;
            
            if 1 % TFE by DCT decomposition                      
               %CutFrq1=[15 10 4.8 2.4 1.2 0.6 0.3 0.1]; 
               %CutFrq1= linspace(0,1,10)*Fs/2;
               
               if 0
                    Npart=10; % equal parts
                    inc=(Fs/2)/Npart;
                    for i=1:Npart-1
                        CutFrq1(i)=Fs/2-i*inc; 
                    end
                end
                if 0 % dyadic filter bank
                    %CutFrq1=[Fs/4 Fs/8 Fs/16 Fs/32 Fs/64];                
                    for i=1:7
                        %CutFrq1(i)=(Fs/2^(i+1)); % dyadic
                        CutFrq1(i)=(Fs/2.1^(i+1)); % dyadic
                    end
                end
             if 0  
                fH(1)=Fs/2; m=1.5; % for dyadic band m=1.5
                for i=1:7
                    fc(i)= ((2*m-1)/(2*m+1))*fH(i);
                    fH(i+1)=fc(i);
                end
                CutFrq1=fc;
                tt=1;               
             end             
             %CutFrq1=0;CutFrq1=[20 15 8.28 6.38 3.967 2 1.0 0.5 0.25]; % based on dips in FFT spectrum
             %CutFrq1(end+1)=0.01;
             %CutFrq1=[.4 .32 .16 0.08 .04 0.02 0.01];
             %CutFrq1=(4:-1:1)*0.1;
               ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(x,t,Fs,CutFrq1,'dct'); 
               figure; 
               for i=1:length(CutFrq1)+1 
                   Data(:,i)=ImfByFilter{i}; 
                   %subplot(length(CutFrq1)+1,1,i);
                   if 1
                   subplot(5,2,i); 
                   plot(t,Data(:,i)); 
                   axis([min(t) max(t) min(Data(:,i)) max(Data(:,i))]);
                   if(i<9)
                       set(gca,'xticklabel',{[]}) 
                   end
                   tmp=int2str(10-i);
                   tmp1=strcat('FIBF',tmp);
                   ylabel(tmp1); FntSize=18;
                   set(gca,'FontSize',FntSize,'FontName','Times')
                 h1 = get(gca, 'xlabel');
                 set(h1,'FontSize',FntSize,'FontName','Times')
                 h1 = get(gca, 'ylabel');
                 set(h1,'FontSize',FntSize,'FontName','Times')
                 end
               end
           end
        if 1 
            figure    
            %plot(sum(Data,2))
            plot(t,x)
            hold on;
            %plot(x)
            plot(t,Data(:,end),'r-.')
            axis([min(t) max(t) min(x) max(x)]);
                   
            xlabel('Year')
            ylabel('Kelvin')
            set(gca,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'xlabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            h1 = get(gca, 'ylabel');
            set(h1,'FontSize',FntSize,'FontName','Times')
            legend('The annual mean global surface temperature anomaly data','Obtained trend of data using proposed FDM');     
        end        
           
        fmax=Fs/2; 
        [freq1 amp1]=sp_DCT_IF_timefrequency_plot(Data,t,Fs,fmax); 
        tt=1;  
                    
    end 
    
    
    
    if 0 % ECG data PLI and BLW removal
        Fs=360;
        x= load('100m.mat');
        t=(0:length(x.val(1,:))-1)/Fs;
        x1=x.val(1,:);
        x2=x.val(2,:);
        subplot(2,1,1); 
        plot(t,x1);
        subplot(2,1,2);       
        plot(t,x2);
        
        range=1:7200;
        y1=x1(range); y1=y1-mean(y1);
        amp=200;
        bl_noise=amp*sin(2*pi*50*t(range))+amp*sin(2*pi*0.1*t(range))+ amp*sin(2*pi*0.3*t(range)) + amp*sin(2*pi*0.45*t(range));
        SNR=10*log10(sum(y1.*y1)/sum(bl_noise.*bl_noise));
        noisysignl=y1+bl_noise;
        subplot(2,1,1); 
        plot(t(range), y1);
        FntSize=18;    
        %set(gca,'xticklabel',{[]})
        set(gca,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'xlabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'ylabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
        %axis([min(t) max(t) min(x) max(x)]);
        subplot(2,1,2);       
        plot(t(range),noisysignl);
        xlabel('Time (s)')
        set(gca,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'xlabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'ylabel');
        set(h1,'FontSize',FntSize,'FontName','Times')        
        tt=1;
        
       CutFrq1=[53 47 0.75];
       %snr=5*p; x=awgn(x,snr,'measured');
       ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(noisysignl,t(range),Fs,CutFrq1,'dct');
       figure;
       for i=1:length(CutFrq1)+1 
           Data(:,i)=ImfByFilter{i}; 
       end
       for i=1:length(CutFrq1)+1      
           
           if 1
               subplot(3,1,1);
               plot(t(range),Data(:,1)+Data(:,3));               
               FntSize=18;
               set(gca,'FontSize',FntSize,'FontName','Times')
               h1 = get(gca, 'xlabel');
               set(h1,'FontSize',FntSize,'FontName','Times')
               h1 = get(gca, 'ylabel');
               set(h1,'FontSize',FntSize,'FontName','Times')
               subplot(3,1,2);
               plot(t(range),Data(:,4));               
               FntSize=18;
               set(gca,'FontSize',FntSize,'FontName','Times')
               h1 = get(gca, 'xlabel');
               set(h1,'FontSize',FntSize,'FontName','Times')
               h1 = get(gca, 'ylabel');
               set(h1,'FontSize',FntSize,'FontName','Times')
               subplot(3,1,3);
               plot(t(range),Data(:,2)); 
               xlabel('Time (s)')
               FntSize=18;
               set(gca,'FontSize',FntSize,'FontName','Times')
               h1 = get(gca, 'xlabel');
               set(h1,'FontSize',FntSize,'FontName','Times')
               h1 = get(gca, 'ylabel');
               set(h1,'FontSize',FntSize,'FontName','Times')
           end
       end
       
    end
    tt=1;
    
    if 1 %1 example from ieee paper
                Fs=8000;  
                t = (0:(1/Fs):1-(1/Fs))';
                x1=(1+0*t.^2).*chirp(t,1000,1,2000)*1;
                %sp_FMD_Low2High_High2LowSacnningGroupDelay(x1,Fs,t); 
                figure
                spectrogram(x1,256,250,256,Fs,'yaxis') ;

                %Fs = 8000; % Sampling rate of signal
                Fc = 780; % Carrier frequency
                %t = [0:Fs-1]'/Fs; % Sampling times
                x2 = 1*sin(2*pi*4*t); % Channel 1
                figure
                plot(x2)
                dev = 200; % Frequency deviation in modulated signal
                y1 = 1*fmmod(x2,Fc,Fs,dev); % Modulate both channels.
                figure
                plot(y1);
                axis([0 length(t) -2 2])
                %z = fmdemod(y,Fc,Fs,dev); % Demodulate both channels.
                %plot(z);
                x=x1+y1;

                figure
                spectrogram(x,256,250,256,Fs,'yaxis') ;
                x=(x'); 
                fmax=Fs/4;
                [freq1 amp1]=sp_DCT_IF_timefrequency_plot(x,t,Fs,fmax);
                figure
                subplot(2,1,1) 
                plot(t,amp1);
                subplot(2,1,2)
                plot(t,freq1);                          
                tt=1; 
            end
    
    if 0 % delta function analysis
        N=1000;
        x=zeros(1,N); 
        x(N/2)=1;
        Fs=100;        
        t=(0:length(x)-1)/Fs;
        fmax=Fs/2;
        [freq1 amp1]=sp_DCT_IF_timefrequency_plot(x,t,Fs,fmax);
        figure
        subplot(2,1,1) 
        plot(t,amp1);
        subplot(2,1,2)
        plot(t,freq1); 
        tt=1;             
    end
    
    [x2,Fs] = audioread('C:\Users\pushpendra.singh\Google Drive\DataBase\speech_F0\cmu_us_bdl_arctic\orig\arctic_a0001.wav');
 
    %dlen=5000:22000;
    dlen=5000:13000;
    t=(0:length(dlen)-1)/Fs;
    t=t';
    x=x2(dlen,1);
    egg=x2(dlen,2);  
       
    L=length(x);   
    t=(0:L-1)/Fs; 
    
    egg=x2(dlen,2);  
    degg=diff(egg); degg=[degg;degg(end)];
    subplot(3,1,1);
    plot(t',x); FntSize=18;
    ylabel('Amplidude')
    %set(gca,'xticklabel',{[]})
    set(gca,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'xlabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'ylabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
    axis([min(t) max(t) min(x) max(x)]);
    subplot(3,1,2);    
    plot(t,egg)
    ylabel('Amplidude')
    %set(gca,'xticklabel',{[]})
    set(gca,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'xlabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'ylabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         
    axis([min(t) max(t) min(egg) max(egg)]);
    subplot(3,1,3);
    plot(t',degg);
    axis([min(t) max(t) min(degg) max(degg)]); FntSize=18;
    xlabel('Time (s)')
    ylabel('Amplidude')
         set(gca,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'xlabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'ylabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
     
    tt=1;
    
    figure;
    %plot(t,x)       
    
    %NFFT = 2^nextpow2(L); % Next power of 2 from length of y 
    NFFT=L;
    Y = fft(x,NFFT)/L;           
    f = Fs/2*linspace(0,1,NFFT/2+1);  
    ampf=2*abs(Y(1:(NFFT/2)+1)); 
    subplot(2,1,1);
    %SigLen=length(f();
    
    [pks,locs] = findpeaks(ampf,'MinPeakDistance',30);
    pkfrq=f(locs);
    plot(f,ampf,f(locs),pks,'or')
    
    plot(f(1:end/3),ampf(1:end/3));  
    %plot(amp(1:130),f(1:130));
    xlabel('f (Hz)')
    ylabel('Speech Spectrum'); FntSize=18;
    set(gca,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'xlabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'ylabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         axis([1 max(f(1:end/3)) min(ampf(1:end/3)) max(ampf(1:end/3))]);
        
    
    Y = fft(degg,NFFT)/L;           
    f = Fs/2*linspace(0,1,NFFT/2+1);  
    ampf=2*abs(Y(1:(NFFT/2)+1)); 
    subplot(2,1,2);
    plot(f(1:end/3),ampf(1:end/3));
    xlabel('f (Hz)')
    ylabel('DEGG Spectrum') 
    set(gca,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'xlabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'ylabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         axis([0 max(f(1:end/3)) min(ampf(1:end/3)) max(ampf(1:end/3))]);
    [val,ind]=max(ampf);
    F0=f(ind); 
    
    TotalNbOfHormonics=12;         
    for i=1:TotalNbOfHormonics
        TmpFrq1(i)=(TotalNbOfHormonics-i)*F0;
    end
    for i=1:TotalNbOfHormonics-1
        TmpFrq2(i)=(TmpFrq1(i)+TmpFrq1(i+1))/2;
    end
    
    %tmp=[flip(pkfrq(1:9)) 0];
    %for i=1:TotalNbOfHormonics-1
    %    TmpFrq2(i)=(tmp(i)+tmp(i+1))/2;
    %end
    
    %for p=1:10  
    if 1 % TFE by DCT decomposition                     
       %CutFrq1=[975 825 675 525 375 225 75];
       CutFrq1=TmpFrq2;
       %snr=5*p; x=awgn(x,snr,'measured');
       ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(x,t,Fs,CutFrq1,'dct');
       figure;
       for i=1:length(CutFrq1)+1 
           Data(:,i)=ImfByFilter{i}; 
           %subplot(length(CutFrq1)+1,1,i);
           if 1
           subplot(6,2,i);
           plot(t,Data(:,i));
           axis([min(t) max(t) min(Data(:,i)) max(Data(:,i))]);
           if(i<11)
               set(gca,'xticklabel',{[]}) 
           end
           tmp=int2str(12-i);
           tmp1=strcat('FIBF',tmp);
           ylabel(tmp1); FntSize=18;
           set(gca,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'xlabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
         h1 = get(gca, 'ylabel');
         set(h1,'FontSize',FntSize,'FontName','Times')
           end
       end
    end
    fmax=4000;
    [freq1 amp1]=sp_DCT_IF_timefrequency_plot(Data,t,Fs,fmax); 
    tt=1; 
    