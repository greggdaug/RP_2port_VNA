% function y=rp_gendut1p(fstart,fstop)


%% Calcualte arbitrary waveform with 16384 samples
% Values of arbitrary waveform must be in range from -1 to 1.

fstart=3e6;
fstop=30e6;
N=16383;
fs=125e6;
ts=1/fs;
T=ts*N;
fbin=fs/N;
freq=0:fs/N:fs/2;
index_fstart=round(fstart/fbin);
index_fstop=round(fstop/fbin);
k=(fstop-fstart)/T;
t=0:ts:N*ts;
x1=sin(2*pi*(fstart.*t+k./2*t.^2));

%% get error correction terms
load('calerrdata1p.mat');

resp=input('Run Continuous? 1=loop, 0=once :');

n=1;

switch resp
    case 1 % continuous
        
        while n==1
    
            %% get dut data
            [signal_num_1,signal_num_2]=rp_getdata(x1);

            a_dut_meas=signal_num_1;
            b_dut_meas=signal_num_2;

            ca_dut_meas=fft(a_dut_meas);
            cb_dut_meas=fft(-b_dut_meas);
            rho_dut_meas=ca_dut_meas./cb_dut_meas;

            for i=index_fstart:index_fstop
                rho_dut_corr(i)=(rho_dut_meas(i)-Ed(i))/((rho_dut_meas(i)*Ep(i)-dE(i)));
                rho_dut_corr_db(i)=20*log10(abs(rho_dut_corr(i)));
                z_dut_corr(i)=50*(1+rho_dut_corr(i))/(1-rho_dut_corr(i));
            end

            figure(1)

            subplot(2,2,1)
            [lineseries,hsm]=smithchart(rho_dut_corr(index_fstart:index_fstop));
            title('Corrected S11');

            subplot(2,2,2)
            plot(freq(index_fstart:index_fstop)/1e6,rho_dut_corr_db(index_fstart:index_fstop),'b')
            title('Corrected Return Loss (dB)');
            grid on

            subplot(2,2,3)
            plot(freq(index_fstart:index_fstop)/1e6,real(z_dut_corr(index_fstart:index_fstop)),'b')
            title('real(Zcorr)');
            grid on

            subplot(2,2,4)
            plot(freq(index_fstart:index_fstop)/1e6,imag(z_dut_corr(index_fstart:index_fstop)),'b')
            title('imag(Zcorr)');
            grid on

%             rho_dut_meas(index_fstart)
% 
%             rho_dut_corr(index_fstart)

        end
        
    case 0 % single
        
        % get dut data
        [signal_num_1,signal_num_2]=rp_getdata(x1);

        a_dut_meas=signal_num_1;
        b_dut_meas=signal_num_2;

        ca_dut_meas=fft(a_dut_meas);
        cb_dut_meas=fft(-b_dut_meas);
        rho_dut_meas=ca_dut_meas./cb_dut_meas;

        for i=index_fstart:index_fstop
            rho_dut_corr(i)=(rho_dut_meas(i)-Ed(i))/((rho_dut_meas(i)*Ep(i)-dE(i)));
            rho_dut_corr_db(i)=20*log10(abs(rho_dut_corr(i)));
            z_dut_corr(i)=50*(1+rho_dut_corr(i))/(1-rho_dut_corr(i));
        end

        figure(1)

        subplot(2,2,1)
        [lineseries,hsm]=smithchart(rho_dut_corr(index_fstart:index_fstop));
        title('Corrected S11');

        subplot(2,2,2)
        plot(freq(index_fstart:index_fstop)/1e6,rho_dut_corr_db(index_fstart:index_fstop),'b')
        title('Corrected Return Loss (dB)');
        grid on

        subplot(2,2,3)
        plot(freq(index_fstart:index_fstop)/1e6,real(z_dut_corr(index_fstart:index_fstop)),'b')
        title('real(Zcorr)');
        grid on

        subplot(2,2,4)
        plot(freq(index_fstart:index_fstop)/1e6,imag(z_dut_corr(index_fstart:index_fstop)),'b')
        title('imag(Zcorr)');
        grid on

%         rho_dut_meas(index_fstart)
% 
%         rho_dut_corr(index_fstart)
   
end

fd=fopen('rpdut.s1p','w');

fprintf(fd,'! File generated with RP 1-port VNA \n');
fprintf(fd,'%s \n','# HZ S RI R 50.0');

for i=index_fstart:index_fstop
    fprintf(fd,'%.0f %.3f %.3f \n',freq(i), real(rho_dut_corr(i)), imag(rho_dut_corr(i)));
end

fclose(fd);

% end
