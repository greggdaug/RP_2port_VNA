% function y=rp_gendut2p(fstart,fstop)


%% Calcualte arbitrary waveform with 16384 samples
% Values of arbitrary waveform must be in range from -1 to 1.

fstart=1e6;
fstop=50e6;
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
load('calerrdata2p.mat');

reply=input('Connect DUT between P1 and P2, press any key to continue');

% resp=input('Run Continuous? 1=loop, 0=once :');

% n=1;

% switch resp
%     case 0 % once

            %% S12 (SW3,SW4 -> COUP2, SW5 -> COUP1)
            
%             pause(0.5);
            
            vna_digio_set(1,1,0);
            
            disp('Calculating S12');

            [signal_num_1,signal_num_2]=rp_getdata(x1);

            s12a_meas=signal_num_1;
            s12b_meas=signal_num_2;

            cs12a_meas=fft(s12a_meas);
            cs12b_meas=fft(s12b_meas);
            rho_s12_meas=cs12a_meas./cs12b_meas;
    
            %% S11 (SW3,SW4,SW5 -> COUP1)
            
            vna_digio_set(0,0,0);
            
            disp('Calculating S11');

            [signal_num_1,signal_num_2]=rp_getdata(x1);

            s11a_meas=signal_num_1;
            s11b_meas=signal_num_2;

            cs11a_meas=fft(s11a_meas);
            cs11b_meas=fft(-s11b_meas);
            rho_s11_meas=cs11a_meas./cs11b_meas;

            %% S22 (SW3,SW4,SW5 -> COUP2)
            
%             pause(0.5);
            
            vna_digio_set(1,1,1);
            
            disp('Calculating S22');

            [signal_num_1,signal_num_2]=rp_getdata(x1);

            s22a_meas=signal_num_1;
            s22b_meas=signal_num_2;

            cs22a_meas=fft(s22a_meas);
            cs22b_meas=fft(-s22b_meas);
            rho_s22_meas=cs22a_meas./cs22b_meas;

            %% S21 (SW3,SW4 -> COUP1, SW5 -> COUP2)
            
%             pause(0.5);
            
            vna_digio_set(0,0,1);
            
            disp('Calculating S21');

            [signal_num_1,signal_num_2]=rp_getdata(x1);

            s21a_meas=signal_num_1;
            s21b_meas=signal_num_2;

            cs21a_meas=fft(s21a_meas);
            cs21b_meas=fft(s21b_meas);
            rho_s21_meas=cs21a_meas./cs21b_meas;



            %% Calculate corrected spar
            
%             for i=index_fstart:index_fstop
% 
%                 D(i)=(1+((rho_s11_meas(i)-e00(i))./e10e01(i)).*e11(i)).*(1+((rho_s22_meas(i)-ep33(i))./ep23ep32(i)).*ep22(i))-((rho_s21_meas(i)-e30(i))./e10e32(i)).*((rho_s12_meas(i)-ep03(i))./ep23ep01(i)).*e22(i).*ep11(i);
%                 
%                 s11_corr(i)=(((rho_s11_meas(i)-e00(i))./e10e01(i)).*(1+((rho_s22_meas(i)-ep33(i))./ep23ep32(i)).*ep22(i))-(e22(i).*((rho_s21_meas(i)-e30(i))./e10e32(i)).*((rho_s12_meas(i)-ep03(i))./ep23ep01(i))))./D(i);
% 
%                 s21_corr(i)=(((rho_s21_meas(i)-e30(i))./e10e32(i)).*(1+((rho_s22_meas(i)-ep33(i))./ep23ep32(i)).*(ep22(i)-e22(i))))./D(i);
% 
%                 s22_corr(i)=(((rho_s22_meas(i)-ep33(i))./ep23ep32(i)).*(1+((rho_s11_meas(i)-e00(i))./e10e01(i)).*e11(i))-(ep11(i).*((rho_s21_meas(i)-e30(i))./e10e32(i)).*((rho_s21_meas(i)-ep03(i))./ep23ep01(i))))./D(i);
% 
%                 s12_corr(i)=(((rho_s12_meas(i)-ep03(i))./ep23ep01(i)).*(1+((rho_s11_meas(i)-e00(i))./e10e01(i)).*(e11(i)-ep11(i))))./D(i);
%             end


%             s11_corr_db=20*log10(abs(s11_corr));
%             s21_corr_db=20*log10(abs(s21_corr));
%             s22_corr_db=20*log10(abs(s22_corr));
%             s12_corr_db=20*log10(abs(s12_corr));

            s11_meas=rho_s11_meas;
            s21_meas=rho_s21_meas;
            s22_meas=rho_s22_meas;
            s12_meas=rho_s12_meas;
            
            s11_meas_db=20*log10(abs(rho_s11_meas));
            s21_meas_db=20*log10(abs(rho_s21_meas));
            s22_meas_db=20*log10(abs(rho_s22_meas));
            s12_meas_db=20*log10(abs(rho_s12_meas));

            figure(1)

%             subplot(2,2,1)
%             [lineseries,hsm]=smithchart(s11_corr(index_fstart:index_fstop));
%             title('Corrected S11');
% 
%             subplot(2,2,2)
%             plot(freq(index_fstart:index_fstop)/1e6,s21_corr_db(index_fstart:index_fstop),'b')
%             title('Corrected S21 (dB)');
%             grid on
% 
%             subplot(2,2,3)
%             plot(freq(index_fstart:index_fstop)/1e6,s12_corr_db(index_fstart:index_fstop),'b')
%             title('Corrected S12 (dB)');
%             grid on
% 
%             subplot(2,2,4)
%             [lineseries,hsm]=smithchart(s22_corr(index_fstart:index_fstop));
%             title('Corrected S22');

            subplot(2,2,1)
            [lineseries,hsm]=smithchart(rho_s11_meas(index_fstart:index_fstop));
            title('Meas S11');

            subplot(2,2,2)
            plot(freq(index_fstart:index_fstop)/1e6,s21_meas_db(index_fstart:index_fstop),'b')
            title('Meas S21 (dB)');
            grid on

            subplot(2,2,3)
            plot(freq(index_fstart:index_fstop)/1e6,s12_meas_db(index_fstart:index_fstop),'b')
            title('Meas S12 (dB)');
            grid on

            subplot(2,2,4)
            [lineseries,hsm]=smithchart(rho_s22_meas(index_fstart:index_fstop));
            title('Meas S22');

        
    
%% Write s2p file
fd=fopen('rpdut2p.s2p','w');

fprintf(fd,'! File generated with RP 2-port VNA \n');
fprintf(fd,'%s \n','# HZ S RI R 50.0');
fprintf(fd,'%s \n','! reS11 imS11 reS21 imS21 reS12 imS12 reS22 imS22');
for i=index_fstart:index_fstop
    fprintf(fd,'%.0f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n',freq(i), real(s11_meas(i)), imag(s11_meas(i)), real(s21_meas(i)), imag(s21_meas(i)), real(s12_meas(i)), imag(s12_meas(i)), real(s22_meas(i)), imag(s22_meas(i)));
    %     fprintf(fd,'%.0f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n',freq(i), real(s11_corr(i)), imag(s11_corr(i)), real(s21_corr(i)), imag(s21_corr(i)), real(s12_corr(i)), imag(s12_corr(i)), real(s22_corr(i)), imag(s22_corr(i)));
end

fclose(fd);

% end
