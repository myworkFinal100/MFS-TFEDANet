function tt=sp_GSTA_DataAnalysis()
close all; 
    clear all;
    format shortEng;
    clc

if 1 %2 example
        if 1
            load gsta.dat;  
            %load LOD_73_14.TXT;
            x=gsta(:,2); 
            t=gsta(:,1);
        end
        if 0
            load gsta1880_2017.txt
            x=gsta1880_2017(:,2); 
            t=gsta1880_2017(:,1);
        end
            
            Fs=1/(t(2)-t(1)); 
            plot(t,x,'LineWidth',2);
            axis([min(t) max(t) min(x) max(x)]);
            grid;  
            title('The annual mean global surface temperature anomaly')
            xlabel('Year')
            ylabel('Kelvin')
            
            
          if 0 % FFT based
            L=length(x);NFFT=2^nextpow2(L);
            Y = fft(x,NFFT)/L; 
            %Y=dct(x,NFFT);
            f = Fs/2*linspace(0,1,NFFT/2+1);  
            ampf=2*abs(Y(1:(NFFT/2)+1)); 
            subplot(2,1,1);
            stem(f,ampf);
            timeScale=1./f(2:end); timeScale=[500 timeScale];
            subplot(2,1,2);
            stem(timeScale,ampf);
            tt=1;
          end
          if 0 % DCT based
            L=length(x);NFFT=2^nextpow2(L);
            %Y=dct(x,NFFT); 
            %NFFT=L;
            Y=dct(x,L);
            f = Fs/2*linspace(0,1,NFFT);  
            ampf=abs(Y(1:(NFFT))); 
            %subplot(2,1,1);
            stem(f,ampf);
            timeScale=1./f(2:end); timeScale=[600 timeScale];
            %subplot(2,1,2);
            %stem((timeScale),(ampf));
            tt=1;   
          end 
           if 1 % TFE by DCT decomposition                     
               %CutFrq1=[.4 .3 .25 .17 0.09 .05]; 
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
                        CutFrq1(i)=(Fs/2^(i+1));
                    end
                end
               
                %CutFrq1(end+1)=0.01;
                %CutFrq1=[.4 .32 .16 0.08 .04 0.02 0.01];
                %CutFrq1=(4:-1:1)*0.1;
                
                
               figure;
               %hold on;
        if 1
               plot(t,x,'LineWidth',1)
               axis([min(t) max(t) min(x) max(x)])
               hold on
               TimeSacle=[35 53 54]; CutFrq1=0; CutFrq=1./TimeSacle; 
               %TimeSacle=[75 100 125]; CutFrq1=0; CutFrq=1./TimeSacle;
               %TimeSacle=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150]; CutFrq1=0; CutFrq=1./TimeSacle;
               %TimeSacle=[10 20 30 40 50 60 70 80]; CutFrq1=0; CutFrq=1./TimeSacle;
               %TimeSacle=[15 25 53 75 100 125]; CutFrq1=0; CutFrq=1./TimeSacle;
               
               for j=1:length(CutFrq)                     
               CutFrq1= [CutFrq(j)];               
               ImfByFilter=0; ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(x,t,Fs,CutFrq1,'dct'); 
               %for i=1:length(CutFrq1)+1 
                   Data=ImfByFilter{2};                    
                   plot(t,Data,'LineWidth',1);                
                       FntSize=18;
                       set(gca,'FontSize',FntSize,'FontName','Times')
                        h1 = get(gca, 'xlabel');
                        set(h1,'FontSize',FntSize,'FontName','Times')
                        h1 = get(gca, 'ylabel');
                        set(h1,'FontSize',FntSize,'FontName','Times')
               end
               legend('GSTA data', '35 years or longer time scale trend','53 years or longer time scale trend','54 years or longer time scale trend');
               
        end      
               
        if 1
               
               TimeSacle=[35 53 54]; CutFrq1=0; CutFrq=1./TimeSacle; 
               %TimeSacle=[75 100 125]; CutFrq1=0; CutFrq=1./TimeSacle;
               %TimeSacle=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150]; CutFrq1=0; CutFrq=1./TimeSacle;
               %TimeSacle=[10 20 30 40 50 60 70 80]; CutFrq1=0; CutFrq=1./TimeSacle;
               %TimeSacle=[15 25 53 75 100 125]; CutFrq1=0; CutFrq=1./TimeSacle;
               
               for j=1:length(CutFrq)                    
               CutFrq1= [CutFrq(j)];               
               ImfByFilter=0; ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(x,t,Fs,CutFrq1,'dct'); 
               %for i=1:length(CutFrq1)+1 
                   Data=ImfByFilter{2}; 
                   %subplot(length(CutFrq1)+1,1,i);
                    variability=x-Data;
                       subplot(3,1,j);
                       if 1
                           plot(t,variability,'LineWidth',1)
                           axis([min(t) max(t) min(variability) max(variability)])
                           tmp=int2str(TimeSacle(j));
                            %tmp1=strcat(tmp, '-years or longer time scale trend');
                            tmp1=strcat(tmp, '-years or longer time scale variability');
                       else
                           plot(t,x,'LineWidth',1)
                           axis([min(t) max(t) min(x) max(x)])
                           hold on 
                           plot(t,Data,'r-.','LineWidth',1); 
                           tmp=int2str(TimeSacle(j));
                           tmp1=strcat(tmp, '-years or longer time scale trend');
                           %tmp1=strcat(tmp, '-years or longer time scale variability');
                       end
                       
                       %plot(t,Data-0.3*(j-1),'LineWidth',2);
                       %hold on
                       %axis([min(t) max(t) min(Data) max(Data)]);
                       if(j<3)
                           set(gca,'xticklabel',{[]}) 
                       end
                       %tmp=int2str(TimeSacle(j));
                       %tmp1=strcat(tmp, '-years or longer time scale trend');
                       %tmp1=strcat(tmp, '-years or longer time scale variability');
                       %legend('GSTA data',tmp1);
                       %ylabel(tmp1);
                       title(tmp1,'FontName','Times');
                       FntSize=18;
                       set(gca,'FontSize',FntSize,'FontName','Times')
                        h1 = get(gca, 'xlabel');
                        set(h1,'FontSize',FntSize,'FontName','Times')
                        h1 = get(gca, 'ylabel');
                        set(h1,'FontSize',FntSize,'FontName','Times')
               end
               %legend('GSTA data', '10 years or longer time scale trend','25 years or longer time scale trend','50 years or longer time scale trend');
               %legend('GSTA data', '75 years or longer time scale trend','100 years or longer time scale trend','125 years or longer time scale trend');
                                    
                   
               %end
               hold off  
               set(gca,'yTickLabel',[]); axis([1850 2010 -3.4 0.5]); 
        end
                
               %TimeSacle=[3 4 6 8 12 16 24 32 44 64 80]; CutFrq1=0; CutFrq1=1./TimeSacle; 
               TimeSacle=[15 25 53 75 100 125]; CutFrq1=0; CutFrq=1./TimeSacle;
               ImfByFilter=sp_DFTOrthogonalOrFIR_IIR_LINOEP(x,t,Fs,CutFrq1,'dct');  
               figure; 
               for i=1:length(CutFrq1)+1 
                   Data(:,i)=ImfByFilter{i}; 
                   %subplot(length(CutFrq1)+1,1,i);
                   if 1
                   subplot(3,2,i); 
                   plot(t,Data(:,i));
                   axis([min(t) max(t) min(Data(:,i)) max(Data(:,i))]);
                   if(i<5)
                       set(gca,'xticklabel',{[]}) 
                   end
                   tmp=int2str(6-i);
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
