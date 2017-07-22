% function y=rp_gencal1p(fstart,fstop)

% 1-port Model
%             
%                     P1      a1
%  a0 o--->---o--->---o--->---o--------
%             |   1   |       |        |
%            e00     e11     S11  DUT  |
%             | e10e01|       |        |
%  b0 o---<---o---<---o---<---o--------
%                             b1
% 
%  e00 = fwd directivity
%  e11 = port1 match
%  (e10e01) = reflection tracking


%% Calcualte arbitrary waveform with 16384 samples
% Values of arbitrary waveform must be in range from -1 to 1.
fstart=3e6;
fstop=30e6;
N=16383;
fs=125e6;
ts=1/fs;
T=ts*N;
fbin=fs/N;
index_fstart=round(fstart/fbin);
index_fstop=round(fstop/fbin);
k=(fstop-fstart)/T;
t=0:ts:N*ts;
x1=sin(2*pi*(fstart.*t+k./2*t.^2));

vna_digio_set(0,0,0);

%% connect short

reply=input('Connect short to P1, press any key to continue');

[signal_num_1,signal_num_2]=rp_getdata(x1);

a_short_meas=signal_num_1;
b_short_meas=signal_num_2;
ca_short_meas=fft(a_short_meas);
cb_short_meas=fft(-b_short_meas);
rho_short_meas=ca_short_meas./cb_short_meas;

% keyboard;
%% connect open
reply=input('Connect open to P1, press any key to continue');

[signal_num_1,signal_num_2]=rp_getdata(x1);

a_open_meas=signal_num_1;
b_open_meas=signal_num_2;
ca_open_meas=fft(a_open_meas);
cb_open_meas=fft(-b_open_meas);
rho_open_meas=ca_open_meas./cb_open_meas;

% keyboard;

%% connect load
reply=input('Connect load to P1, press any key to continue');

[signal_num_1,signal_num_2]=rp_getdata(x1);

a_load_meas=signal_num_1;
b_load_meas=signal_num_2;
ca_load_meas=fft(a_load_meas);
cb_load_meas=fft(-b_load_meas);
rho_load_meas=ca_load_meas./cb_load_meas;

%% SOL factors
% rho_load=complex(-0.005*ones(1,N+1),0.1*ones(1,N+1));
rho_load=complex(-0.001*ones(1,N+1),0.001*ones(1,N+1));
rho_short=complex(-1*ones(1,N+1),0.09*ones(1,N+1));
rho_open=complex(0.98*ones(1,N+1),0.08*ones(1,N+1));
% rho_open=complex(-1*ones(1,N+1),-1*ones(1,N+1));
eones=complex(ones(1,N+1),ones(1,N+1));

%% generate error terms
for i=index_fstart:index_fstop
    errmat(:,i)=inv([[eones(i),rho_short(i)*rho_short_meas(i),-rho_short(i)];[eones(i),rho_open(i)*rho_open_meas(i),-rho_open(i)];[eones(i),rho_load(i)*rho_load_meas(i),-rho_load(i)]])*[rho_short_meas(i);rho_open_meas(i);rho_load_meas(i)];
    Ed(i)=errmat(1,i); % e00
    Ep(i)=errmat(2,i); % e11
    dE(i)=errmat(3,i); % e10e01
%     rho_dut(i)=rho_dut_meas(i);
%     rho_dut_db(i)=20*log10(abs(rho_dut(i)));
%     rho_dut_corr(i)=(rho_dut(i)-Ed(i))/((rho_dut(i)*Ep(i)-dE(i)));
%     z_dut_corr(i)=50*(1+rho_dut_corr(i))/(1-rho_dut_corr(i));
%     rho_dut_corr_db(i)=20*log10(abs(rho_dut_corr(i)));
end

save('calerrdata1p.mat','Ed','Ep','dE');
save('calmeasdata1p.mat','rho_short_meas','rho_open_meas','rho_load_meas');

% end


