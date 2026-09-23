% References:   
% (1) Breaking the Limits ? Redefining the Instantaneous Frequency, http://arxiv.org/abs/1605.00975
% (2) The Fourier Decomposition Method for nonlinear and nonstationary time
% series analysis, http://arxiv.org/abs/1503.06675
% To obtain FIBFs using FDM based on 'dct' or 'dft' or 'fir' use function sp_DFTOrthogonalOrFIR_IIR_LINOEP

function sp_Main_examples()
    close all; 
    clear all;    
    clc;  
    
    if 0
        load quadchirp;
        fs = 1000;
        [S,F,T] = spectrogram(quadchirp,100,98,128,fs);
        helperCWTTimeFreqPlot(S,T,F,'surf','STFT of Quadratic Chirp','Seconds','Hz')
    end
    
            
    if 1
        if 1 % five parallel chirp signal
            Fs=1*8000; 
            t = 0:(1/Fs):1-(1/Fs);
            t=t';
            x0=(1+0*t.^2).*chirp(t,500,1,1500)*1; 
            x1=(1+0*t.^2).*chirp(t,1000,1,2000)*1;
            x2=(1+0*t.^2).*chirp(t,1500,1,2500)*1;
            x3=(1+0*t.^2).*chirp(t,2000,1,3000)*1;
            x4=(1+0*t.^2).*chirp(t,2500,1,3500)*1;
            x=x0+x1+x2+x3+x4;

            if 0 % two delayed chirps
               x1=[x0; zeros(length(x0)/2,1)]; 
               x2=[zeros(length(x0)/2,1);x0]; 
               t=[t; t(end)+t(1)+t(1:end/2)];
               x=x2+x1;
            end

            Fmax=Fs/2; 
            Npart=20; 
        end 
        if 0 %1 example from ieee paper, sum of FM+chipr
            Fs=8000;  
            t = (0:(1/Fs):1-(1/Fs))';
            x1=chirp(t,1000,1,2000)*1;

            Fmax=Fs/4;                
            [freq amp]=sp_IF_timefrequency_plot(x1,t,Fs,Fmax); % plot without decomposition

            %figure
            %spectrogram(x1,256,250,256,Fs,'yaxis') ;

            Fc = 780; % Carrier frequency               
            x2 = 1*sin(2*pi*4*t); % Channel 1
            figure
            plot(x2)
            dev = 200; % Frequency deviation in modulated signal
            y1 = 1*fmmod(x2,Fc,Fs,dev); % Modulate both channels.

            [freq amp]=sp_IF_timefrequency_plot(y1,t,Fs,Fmax); % plot without decomposition                
            x=x1+y1;
            Npart=100; 

        end
        [freq amp]=sp_IF_timefrequency_plot(x,t,Fs,Fmax); % plot without decomposition
                  
        %CutFrq1=[0 Fs/10 Fs/9 Fs/8 Fs/7 Fs/5 Fs/4 Fs/3 Fs/2]; % Manual desired cutoff frequencies
        [CutFrq1] = sp_SetCutOffFreq(x,Fs,'ub',10); %uniform filter bank 10
        %[CutFrq1] = sp_SetCutOffFreq(x,Fs,'d',10); %dyadic filter bank 10
        %[CutFrq1] = sp_SetCutOffFreq(x,Fs,'elp',10,1); %equal Lp-norm band 10
        
        %To obtain FIBFs using FDM: ('dct', 'dft', 'fir') use function sp_DFTOrthogonalOrFIR_IIR_LINOEP  
        ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(x,t,Fs,CutFrq1,'dct'); % 'dft' or 'fir' 

        for i=1:length(CutFrq1)-1 
            Data(:,i)=ImfByFilter{i}; 
        end                            
        [freq amp]=sp_IF_timefrequency_plot(Data,t,Fs,Fmax);                                       
                  
        tt=1;                                                        
    end 
        
close all;        
if 1 % delta function
    N=4000;
    Fs=1000;
    t=(0:N-1)/Fs;
    t=t;
    x=zeros(1,N);
    %X(N/2)=1;
    x((N/2)-1)=1; 
    plot(0:length(x)-1,x);           
    Fmax=Fs/2;
    %CutFrq1=[Fs/4 Fs/8 Fs/16 Fs/32 Fs/64]; % dydic FB
    [freq amp]=sp_IF_timefrequency_plot(x,t,Fs,Fmax);


    [CutFrq1] = sp_SetCutOffFreq(x,Fs,'ub',10); %uniform filter bank 10
    ImfByFilter=[];
    ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(x,t,Fs,CutFrq1,'dft'); % use 'dft' or 'fir'

    Data=[];
    for i=1:length(CutFrq1)-1
        Data(:,i)=ImfByFilter{i};
    end

    [freq amp]=sp_IF_timefrequency_plot(Data,t,Fs,Fmax); 
    tt=1;            
end