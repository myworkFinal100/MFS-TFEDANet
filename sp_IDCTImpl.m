function x = sp_IDCTImpl(y)
    N = length(y);
    if(N == 1)
        x = y(1);
    else
        Q=exp(pi*1i*((0:(N-1))')/(2*N));
        Q(2:N)=Q(2:N)/sqrt(2);
        yrev=y(N:(-1):2);
        toapply=[ y(1); Q(2:N).*(y(2:N)-1i*yrev) ];
        %x1=IFFTImpl(toapply);
        x1=ifft(toapply)*sqrt(N);
        x=zeros(N,1);
        x(1:2:(N-1))=x1(1:(N/2));
        x(2:2:N)=x1(N:(-1):(N/2+1));
        x=real(x);
    end

    