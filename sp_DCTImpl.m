function y = sp_DCTImpl(x)
    N = length(x);
    if(N == 1)
        y = x;
    else
        x1 = [x(1:2:(N-1)); x(N:(-2):2)];
        %y = FFTImpl(x1);
        y=fft(x1)/sqrt(N);
        rp = real(y);
        ip = imag(y);
        y = cos(pi*((0:(N-1))')/(2*N)).*rp + sin(pi*((0:(N-1))')/(2*N)).*ip;
        y(2:N) = sqrt(2)*y(2:N);
        y=real(y);
    end
    