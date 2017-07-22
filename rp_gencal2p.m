% function y=rp_gencal2p(fstart,fstop)

% RedPitaya 2port VNA Calibration
% 12-term SOL Cal
% Re: Network Analyzer Error Models and Calibration Methods by Doug Rytting (Agilent)
% G. Daugherty 3/4/17
% 
% Forward Model
%                                e30
%              ------------------->-------------------
%             |                                       |
%             |                  DUT                  |
%             |       P1      a1      b2      P2      |
%  a0 o--->---o--->---o--->---o--->---o--->---o--->---o--->---o b3
%             |   1   |       |  S21  |       | e10e32
%            e00     e11     S11     S22     e22
%             | e10e01|       |  S12  |       |
%  b0 o---<---o---<---o---<---o---<---o---<---o
%                             b1      a2
% 
%  e00 = directivity
%  e11 = port1 match
%  (e10e01) = reflection tracking
%  (e10e32) = transmission tracking
%  e22 = port2 match
%  e30 = leakage
% 
% 
%  Reverse Model
% 
%                                DUT
%                     P1     a'1     b'2      P2
%                     o--->---o--->---o--->---o--->---o--->---o b'3
%                     |       |  S21  |       | e'23e'32
%                    e'11    S11     S22     e'22    e'33
%                     |       |  S12  |       |   1   |
%  b'0 o---<---o---<--o---<---o---<---o---<---o---<---o---<---o a'3
%               e'23e'01     b'1     a'2              |
%              |                                      |
%              |                                      |
%               ------------------<-------------------
%                                e'03
% 
%  e'33 = directivity
%  e'11 = port1 match
%  (e'23e'32) = reflection tracking
%  (e'23e'01) = transmission tracking
%  e'22 = port2 match
%  e'03 = leakage
%
%  Switch Info
% 
%    SW3=DIO3
%    SW4=DIO4
%    SW5=DIO5
%                              |-------|
%               |-----|o------>| COUP2 |---------------> P2
%  RP DAC  o--->| SW3 |        |-------|
%               |-----|o         o   o     |-------|
%                      | --------|---|---->| COUP1 |---> P1
%                                |   |     |-------|
%                                |   |       o   o
%                                |   ---|    |   |
%                                |      |    |   |
%                                o   o<-|-----   |
%                               |-----| |        |
%  RP ADC2 o<-------------------| SW4 | |        |
%                               |-----| |        |
%                                       |        |
%                                |-------        |
%                                o   o<-----------
%                               |-----|
%  RP ADC1 o<-------------------| SW5 |
%                               |-----|
%                               
%                              
%
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


%% SOL factors
rho_load=complex(-0.005*ones(1,N+1),0.1*ones(1,N+1));
rho_short=complex(-1*ones(1,N+1),0.09*ones(1,N+1));
rho_open=complex(0.98*ones(1,N+1),0.08*ones(1,N+1));
% rho_open=complex(-1*ones(1,N+1),-1*ones(1,N+1));
eones=complex(ones(1,N+1),ones(1,N+1));


%% P1(Fwd) short (SW3,SW4,SW5 -> COUP1)

vna_digio_set(0,0,0);

reply=input('Connect short to P1, press any key to continue');

disp('Calculating Fwd Short');

[signal_num_1,signal_num_2]=rp_getdata(x1);

p1a_short_meas=signal_num_1;
p1b_short_meas=signal_num_2;
cp1a_short_meas=fft(p1a_short_meas);
cp1b_short_meas=fft(-p1b_short_meas);
rho_meas_short_fwd=cp1a_short_meas./cp1b_short_meas;

[lineseries,hsm]=smithchart(rho_meas_short_fwd(index_fstart:index_fstop));

keyboard;
%% P1(Fwd) open (SW3,SW4,SW5 -> COUP1)

reply=input('Connect open to P1, press any key to continue');

disp('Calculating Fwd Open');

[signal_num_1,signal_num_2]=rp_getdata(x1);

p1a_open_meas=signal_num_1;
p1b_open_meas=signal_num_2;
cp1a_open_meas=fft(p1a_open_meas);
cp1b_open_meas=fft(-p1b_open_meas);
rho_meas_open_fwd=cp1a_open_meas./cp1b_open_meas;

[lineseries,hsm]=smithchart(rho_meas_open_fwd(index_fstart:index_fstop));

keyboard;
%% P1(Fwd) load (SW3,SW4,SW5 -> COUP1)

reply=input('Connect load to P1, press any key to continue');

disp('Calculating Fwd Load');

[signal_num_1,signal_num_2]=rp_getdata(x1);

p1a_load_meas=signal_num_1;
p1b_load_meas=signal_num_2;
cp1a_load_meas=fft(p1a_load_meas);
cp1b_load_meas=fft(-p1b_load_meas);
rho_meas_load_fwd=cp1a_load_meas./cp1b_load_meas;

[lineseries,hsm]=smithchart(rho_meas_load_fwd(index_fstart:index_fstop));

keyboard;
%% Fwd leakage e30 (SW3,SW4 -> COUP1, SW5 -> COUP2)

vna_digio_set(0,0,1);

reply=input('Connect loads to P1 and P2, press any key to continue');

disp('Calculating Fwd Leakage');

[signal_num_1,signal_num_2]=rp_getdata(x1);

e30a_meas=signal_num_1;
e30b_meas=signal_num_2;
ce30a_meas=fft(e30a_meas);
ce30b_meas=fft(-e30b_meas);
e30=ce30a_meas./ce30b_meas;

[lineseries,hsm]=smithchart(e30(index_fstart:index_fstop));

keyboard;
%% Fwd s11m (SW1,SW2,SW3 -> COUP1)

vna_digio_set(0,0,0);

reply=input('Connect thru P1 -> P2, press any key to continue');

disp('Calculating Fwd S11m');

[signal_num_1,signal_num_2]=rp_getdata(x1);

s11ma_meas=signal_num_1;
s11mb_meas=signal_num_2;
cs11ma_meas=fft(s11ma_meas);
cs11mb_meas=fft(-s11mb_meas);
s11m=cs11ma_meas./cs11mb_meas

[lineseries,hsm]=smithchart(s11m(index_fstart:index_fstop));

keyboard;
%% Fwd s21m (SW3,SW4 -> COUP1, SW5 -> COUP2)

vna_digio_set(0,0,1);

reply=input('Connect thru P1 -> P2, press any key to continue');

disp('Calculating Fwd S21m');

[signal_num_1,signal_num_2]=rp_getdata(x1);

s21ma_meas=signal_num_1;
s21mb_meas=signal_num_2;
cs21ma_meas=fft(s21ma_meas);
cs21mb_meas=fft(-s21mb_meas);
s21m=s21ma_meas./s21mb_meas;

[lineseries,hsm]=smithchart(s21m(index_fstart:index_fstop));

keyboard;
%% generate fwd error terms
for i=index_fstart:index_fstop
    errmatfwd(:,i)=inv([[eones(i),rho_short(i)*rho_meas_short_fwd(i),-rho_short(i)];[eones(i),rho_open(i)*rho_meas_open_fwd(i),-rho_open(i)];[eones(i),rho_load(i)*rho_meas_load_fwd(i),-rho_load(i)]])*[rho_meas_short_fwd(i);rho_meas_open_fwd(i);rho_meas_load_fwd(i)];
    e00(i)=errmatfwd(1,i); % e00
    e11(i)=errmatfwd(2,i); % e11
    dEfwd(i)=errmatfwd(3,i); % delta
    s11mf(i)=s11m(i);
    s21mf(i)=s21m(i);
    e30f(i)=e30(i);
end

e10e01=(e00.*e11)-dEfwd;
e22=(s11mf-e00)./(s11mf.*e11-dEfwd);
e10e32=(s21mf-e30f)./(1-e11.*e22);

keyboard;
%% P2(Rev) short (SW3,SW4,SW5 -> COUP2)

vna_digio_set(1,1,1);

reply=input('Connect short to P2, press any key to continue');

disp('Calculating Rev Short');

[signal_num_1,signal_num_2]=rp_getdata(x1);

p2a_short_meas=signal_num_1;
p2b_short_meas=signal_num_2;
cp2a_short_meas=fft(p2a_short_meas);
cp2b_short_meas=fft(-p2b_short_meas);
rho_meas_short_rev=cp2a_short_meas./cp2b_short_meas;

[lineseries,hsm]=smithchart(rho_meas_short_rev(index_fstart:index_fstop));

keyboard;
%% P2(Rev) open (SW3,SW4,SW5 -> COUP2)

reply=input('Connect open to P2, press any key to continue');

disp('Calculating Rev Open');

[signal_num_1,signal_num_2]=rp_getdata(x1);

p2a_open_meas=signal_num_1;
p2b_open_meas=signal_num_2;
cp2a_open_meas=fft(p2a_open_meas);
cp2b_open_meas=fft(-p2b_open_meas);
rho_meas_open_rev=cp2a_open_meas./cp2b_open_meas;

[lineseries,hsm]=smithchart(rho_meas_open_rev(index_fstart:index_fstop));

keyboard;
%% P2(Rev) load (SW3,SW4,SW5 -> COUP2)

reply=input('Connect load to P1, press any key to continue');

disp('Calculating Rev Load');

[signal_num_1,signal_num_2]=rp_getdata(x1);

p2a_load_meas=signal_num_1;
p2b_load_meas=signal_num_2;
cp2a_load_meas=fft(p2a_load_meas);
cp2b_load_meas=fft(-p2b_load_meas);
rho_meas_load_rev=cp2a_load_meas./cp2b_load_meas;

[lineseries,hsm]=smithchart(rho_meas_load_rev(index_fstart:index_fstop));

keyboard;
%% Rev leakage ep03 (S12) (SW3,SW4 -> COUP2 SW5 -> COUP1)

vna_digio_set(1,1,0);

reply=input('Connect loads to P1 -> P2, press any key to continue');

disp('Calculating Rev Leakage');

[signal_num_1,signal_num_2]=rp_getdata(x1);

ep03a_meas=signal_num_1;
ep03b_meas=signal_num_2;
cep03a_meas=fft(ep03a_meas);
cep03b_meas=fft(-ep03b_meas);
ep03=cep03a_meas./cep03b_meas;

[lineseries,hsm]=smithchart(ep03(index_fstart:index_fstop));

keyboard;
%% Rev s22m (SW3,SW4,SW5 -> COUP2)

vna_digio_set(1,1,1);

reply=input('Connect thru P1 -> P2, press any key to continue');

disp('Calculating Rev S22m');

[signal_num_1,signal_num_2]=rp_getdata(x1);

s22ma_meas=signal_num_1;
s22mb_meas=signal_num_2;
cs22ma_meas=fft(s22ma_meas);
cs22mb_meas=fft(-s22mb_meas);
s22m=cs22ma_meas./cs22mb_meas;

[lineseries,hsm]=smithchart(s22m(index_fstart:index_fstop));

keyboard;
%% Rev s12m (SW3,SW4 -> COUP2, SW5 -> COUP1)

vna_digio_set(1,1,0);

disp('Calculating Rev S12m');

reply=input('Connect thru P1 -> P2, press any key to continue');

[signal_num_1,signal_num_2]=rp_getdata(x1);

s12ma_meas=signal_num_1;
s12mb_meas=signal_num_2;
cs12ma_meas=fft(s12ma_meas);
cs12mb_meas=fft(-s12mb_meas);
s12m=cs12ma_meas./cs12mb_meas;

[lineseries,hsm]=smithchart(s12m(index_fstart:index_fstop));

keyboard;
%% generate rev error terms
for i=index_fstart:index_fstop
    errmatrev(:,i)=inv([[eones(i),rho_short(i)*rho_meas_short_rev(i),-rho_short(i)];[eones(i),rho_open(i)*rho_meas_open_rev(i),-rho_open(i)];[eones(i),rho_load(i)*rho_meas_load_rev(i),-rho_load(i)]])*[rho_meas_short_rev(i);rho_meas_open_rev(i);rho_meas_load_rev(i)];
    ep33(i)=errmatrev(1,i); % ep33
    ep22(i)=errmatrev(2,i); % ep22
    dErev(i)=errmatrev(3,i); % delta
    s22mf(i)=s22m(i);
    s12mf(i)=s12m(i);
    ep03f(i)=ep03(i);
end

ep23ep32=(ep33.*ep22)-dErev;
ep11=(s22mf-ep33)./(s22mf.*ep22-dEfwd);
ep23ep01=(s12mf-ep03f).*(1-(ep33.*ep22));

keyboard;
%% save 12 error terms

disp('Saving 2p error terms');

save('calerrdata2p.mat','e00','e11','e22','e30','e10e01','e10e32','ep11','ep22','ep33','ep03','ep23ep01','ep23ep32');
% save('calmeasdata2p.mat','rho_meas_short_fwd','rho_meas_open_fwd','rho_meas_load_fwd');

% end


