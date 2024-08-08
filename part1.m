%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate stream of random bits (10,000 bit)
N = 10000;                  % number of bits
bits = randi([0 1], 1, N);  % random binary vector

%% Pick a type and assign it to the variable "type"
types = ["UPNRZ"; "UPRZ"; "PNRZ"; "BPRZ"; "MANCHESTER"];
signal_type = "BPRZ";

%% Choose a bit rate (bit/s)
Rb = 1000;
Tb = 1/Rb;

%% Sampling frequency
fs = 30000;             % 30KHz
ts = 1/fs;              % sampling period
ds = fs/Rb;             % samples per bit

%% time domain
T = N * Tb;             % total time
t = 0 : ts : (T-ts);    % time vector

flag = 0;               % flag for invalid signal_type input
global se;              % count for BPRZ error detection circuit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p_all(time_vector, amplitude_vector,
              samples_per_bit, bit_period, n)
    % Plot time domain
    figure;
    plot(time_vector, amplitude_vector);
    zoom xon;
    xlabel('Time (s)');
    ylabel('Amplitude(V)');

    % Plot the spectral domains of the pulses
    figure;
    Y = fft(amplitude_vector);      % Fourier transform of the signal
    Y = Y .* conj(Y);               % square of the magnitude
    f = linspace(-0.5, 0.5, n*samples_per_bit);  % frequency vector
    plot(f, fftshift(Y));           % plot frequency spectrum
    xlabel('Normalized frequency');
    ylabel('Power spectrum');

    % Eye diagram
    eyediagram(amplitude_vector,3,bit_period);
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(signal_type, "UPNRZ") == 1
  k = 1;
  for i = 1:N
    if bits(i) == 1
      for j = 1:ds
        y(k) = 1.2;
        k = k + 1;
      endfor
    elseif bits(i) == 0
      for j = 1:ds
        y(k) = 0;
        k = k + 1;
      endfor
    end
  end

elseif strcmp(signal_type, "PNRZ") == 1
  k = 1;
  for i = 1:N
    if bits(i) == 1
      for j = 1:ds
        y(k) = 1.2;
        k = k + 1;
      endfor
    elseif bits(i) == 0
      for j = 1:ds
        y(k) = -1.2;
        k = k + 1;
      endfor
    end
  end

elseif strcmp(signal_type, "UPRZ") == 1
  k = 1;
  for i = 1:N
    if bits(i) == 1
      for j = 1:0.5 * ds
        y(k) = 1.2;
        k = k + 1;
      endfor
      for j = (0.5 * ds)+1 : ds
        y(k) = 0;
        k = k + 1;
      endfor
    elseif bits(i) == 0
      for j = 1:ds
        y(k) = 0;
        k = k + 1;
      endfor
    end
  endfor

elseif  strcmp(signal_type, "BPRZ") == 1
  k = 1;
  sign = 1;
  for i = 1:N
    if bits(i) == 1
      for j = 1:0.5 * ds
        y(k) = sign * 1.2;
        k = k + 1;
      endfor
      sign = sign * (-1);
      for j = (0.5 * ds)+1 : ds
        y(k) = 0;
        k = k + 1;
      endfor
    elseif bits(i) == 0
      for j = 1:ds
        y(k) = 0;
        k = k + 1;
      endfor
    end
  endfor

elseif strcmp(signal_type, "MANCHESTER") == 1
  k = 1;
  for i = 1:N
    if bits(i) == 1
      for j = 1:0.5 * ds
        y(k) = 1.2;
        k = k + 1;
      endfor
      for j = (0.5 * ds)+1 : ds
        y(k) = -1.2;
        k = k + 1;
      endfor
    elseif bits(i) == 0
      for j = 1:0.5 * ds
        y(k) = -1.2;
        k = k + 1;
      endfor
      for j = (0.5 * ds)+1 : ds
        y(k) = 1.2;
        k = k + 1;
      endfor
    end
  endfor

else
  flag = 1;
  printf("Unsupported type of line coding.\n");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Decision Device %%%%%%%%%%%%%%%%%%%%%%%%%%

function rxb = decision_device(n,sig_type,samples_per_bit,
                              amplitude_vector)
  global se;

  if strcmp(sig_type, "UPNRZ") == 1
    k = 1;
    c = 1;
    for i = 1:n
        total = 0;
        for j = 1:samples_per_bit
          total = total + amplitude_vector(k);
          k = k + 1;
        endfor

        avg = total/samples_per_bit;
        if avg < 0.6
          rxb(i) = 0;
        else
          rxb(i) = 1;
        end
    endfor

  elseif strcmp(sig_type, "PNRZ") == 1
    k = 1;
    for i = 1:n
        total = 0;
        for j = 1:samples_per_bit
          total = total + amplitude_vector(k);
          k = k + 1;
        endfor
        avg = total/samples_per_bit;
        if avg < 0
          rxb(i) = 0;
        else
          rxb(i) = 1;
        end
    endfor

  elseif strcmp(sig_type, "UPRZ") == 1
    k = 1;
    for i = 1:n
      total1 = 0;
      total2 = 0;
      for j = 1:0.5 * samples_per_bit
        total1 = total1 + amplitude_vector(k);
        k = k + 1;
      endfor
      avg1 = total1/ (0.5*samples_per_bit);
      for o = (0.5*samples_per_bit) +1 :samples_per_bit
        total2 = total2 + amplitude_vector(k);
        k = k + 1;
      endfor
      avg2 = total2 / (0.5*samples_per_bit);
      if (avg1 > 0.6)
        rxb(i) = 1;
      elseif (avg1 < 0.6)
        rxb(i) = 0;
      endif
    endfor

  elseif strcmp(sig_type, "BPRZ") == 1
    k = 1;
    for i = 1:n
      total1 = 0;
      total2 = 0;
      sign = 1;
      se = 0;
      for j = 1:0.5 * samples_per_bit
        total1 = total1 + amplitude_vector(k);
        k = k + 1;
      endfor
      avg1 = total1/ (0.5*samples_per_bit);
      for o = (0.5*samples_per_bit) +1 :samples_per_bit
        total2 = total2 + amplitude_vector(k);
        k = k + 1;
      endfor
      avg2 = total2 / (0.5*samples_per_bit);
      if (abs(avg1) > 0.6)
        rxb(i) = 1;

%% error detection circuit (bonus) %%
        if avg1*sign > 0.6
          sign = sign*(-1);
        else
          se = se + 1;
        endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      elseif (abs(avg1) < 0.6)
        rxb(i) = 0;
      endif
    endfor

  elseif strcmp(sig_type, "MANCHESTER") == 1
    k = 1;
    for i = 1:n
      total1 = 0;
      total2 = 0;
      for j = 1:0.5 * samples_per_bit
        total1 = total1 + amplitude_vector(k);
        k = k + 1;
      endfor
      avg1 = total1/ (0.5*samples_per_bit);
      for o = (0.5*samples_per_bit) +1 :samples_per_bit
        total2 = total2 + amplitude_vector(k);
        k = k + 1;
      endfor
      avg2 = total2 / (0.5*samples_per_bit);
      if (avg1 > 0.6)
        rxb(i) = 1;
      elseif (avg1 < -0.6)
        rxb(i) = 0;
      endif
    endfor

  end
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = noise(time_vector, amplitude_vector, sigma)
  noise = sigma * randn(1, length(time_vector));
  s = amplitude_vector + noise;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = calculate_BER(n,txb,rxb)
  error_count = 0;
  for i = 1:n
    if rxb(i) != txb(i)
      error_count = error_count + 1;
    endif
  endfor
  val = error_count / n;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sigma Sweep %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BER_s = sigma_sweep(t_vector,a_vector,n,sig_type,
                            samples_per_bit,txb)
  sigma2 = [0.1 : 0.1 : 1];
  for sc = 1:10
    error_count = 0;
    s = noise(t_vector, a_vector, sigma2(sc));
    sweep_bits = decision_device(n,sig_type,samples_per_bit,s);
    BER_s(sc) = calculate_BER(n,txb,sweep_bits);
  endfor

  figure;
  semilogy(BER_s, sigma2);
  xlabel('BER');
  ylabel('Noise RMS');
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag == 0
  %% No noise
  received_bits1 = decision_device(N, signal_type, ds, y);
  BER_no_noise = calculate_BER(N, bits, received_bits1);
  p_all(t, y, ds, Tb, N);

  %% With Noise
  s = noise(t, y, 0.1);
  received_bits2 = decision_device(N, signal_type, ds, s);
  BER_with_noise = calculate_BER(N, bits, received_bits2);
  p_all(t, s, ds, Tb, N);

  %% Sigma Sweep
  BER_sweep = sigma_sweep(t, y, N, signal_type, ds, bits);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% End Of Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
