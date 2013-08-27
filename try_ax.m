%ultrasonic_params;
fs=450450;
% Nfft=[0.001 0.0005 0.00025];
fft_length=0.0005;  % s
NW=22;
K=43;
pval=0.01;

base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_sys_test';
date_str='06132012';
letter_str='D';
data_analysis_dir_name='Data_analysis10';
input_files_base_name=fullfile(base_dir_name, ...
                               ['sys_test_' date_str], ...
                               'demux', ...
                               ['Test_' letter_str '_1']);
output_file_name='output.ax';

t_start=23;  % second that contains voc 51
t_stop=24;

ax1(fs,fft_length,NW,K,pval,input_files_base_name,output_file_name,t_start,t_stop)
