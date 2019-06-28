
Fs = length(Whiskers_time_0)/((Whiskers_time_0(end)-Whiskers_time_0(1))/1000);
Fn = Fs/2;                                              % Nyquist Frequency
Wp = 6/Fn;                                        % Passband Frequencies (Normalized)
Ws = 8/Fn;                                        % Stopband Frequencies (Normalized)
Rp = 10;                                                % Passband Ripple (dB)
Rs = 50;                                                % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Filter Order
[c,b,a] = cheby2(n,Rs,Ws);                              % Filter Design
[sosbp,gbp] = zp2sos(c,b,a);                            % Convert To Second-Order-Section For Stability
figure
freqz(sosbp, 2^16, Fs)                                  % Bode Plot Of Filter



WSP = filtfilt(sosbp, gbp, Whiskers_angle_0);

figure, plot(WSP)
%% Get rid of nan
%
not_nan = find(~isnan(Whiskers_angle_0));
Whiskers_angle_0 = Whiskers_angle_0(not_nan);
Whiskers_time_0 = Whiskers_time_0(not_nan);

%%