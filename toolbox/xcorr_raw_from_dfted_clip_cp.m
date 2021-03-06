function [xcorr_raw,tau_line]=xcorr_raw_from_dfted_clip(V,dt,a,M,verbosity)

% a is the gain for each channel, and these are multiplied in.

% colors for microphones, if needed
clr_mike=[0 0   1  ; ...
          0 0.7 0  ; ...
          1 0   0  ; ...
          0 0.8 0.8];

% calculate the time lag for each element of xcorr_raw
[N,K]=size(V);  % K the number of mikes
tau_line=fftshift(fft_base(N,dt));  % want large neg times first

% calculate the cross power spectrum for each pair, show
n_pairs=K*(K-1)/2;
xcorr_raw=zeros(N,n_pairs);
for i_pair=1:n_pairs
  i_mike_pair_this=find(M(i_pair,:));
  i_mike=i_mike_pair_this(1);
  j_mike=i_mike_pair_this(2);

  % calc cross-power spectrum
  Xcorr_raw_this= ...
    a(i_mike)*a(j_mike)*(V(:,i_mike).*conj(V(:,j_mike)));

  % go to the time domain
  xcorr_raw_this=fftshift(real(ifft(Xcorr_raw_this)));  % want large neg times first

  % plot that
  if verbosity>=3
    T_max=0.002;  % s
    title_str_this= ...
      sprintf('Cross correlation function between mic %d and %d', ...
              i_mike,j_mike);
    figure('color','w');
    plot(1000*tau_line,1e12*xcorr_raw_this, ...
         'color',mean([clr_mike(i_mike,:);clr_mike(j_mike,:)]));
    xlabel('Lag (ms)');
    ylabel('Cross-correlation (mV^4)');
    title(title_str_this,'interpreter','none');
    xlim(1000*[-T_max T_max]);            
    drawnow;
  end

  % store xcorrs
  xcorr_raw(:,i_pair)=xcorr_raw_this;
end  

end
