function [ video_pulse_start_ts dist] = fn_video_pulse_start_ts(dir1, dir2, audio_fname_prefix, video_fname_prefix, load_time_stamps, fc, vfc) 
%fn_video_pulse_start_ts
%   function extracts timestamps associated with recorded pulses from
%   function generator that drive camera frames
%
%   OUTPUT (video_pulse_start_ts) is the timestamps associated with the
%   start of the frame
%
%   will plot voltages and corresponding start of pulses (red dots)
%   will plot distribution of timestamps in recorded steampix video (info
%   from fname seq file)
%
%   Variables 
%
%   dir1 = name of directory with audio data
%   dir2 = name of directory with streampix video file
%   audio_fname_prefix = name of saved audio file
%   video_fname_prefix = name of saved video file
%   load_time_stamps = 'y' or 'n'  and will either load saved data or
%   create file to save

if strcmp(load_time_stamps,'n')==1
%     cd (dir1)
%     filename1 = sprintf('%s.ch5',audio_fname_prefix);
%     fid = fopen(filename1);
%     tmp = fread(fid,'float32');
    
    %modified on 11/2/12 to load in large files--does it one chunck at a
    %time.
%     tmp_a2 = tmp;
%     tmp_i2=tmp_a2>1;
%     video_pulse_start_ts = find(diff(tmp_i2)==1);
    video_pulse_start_ts = fn_video_time_ts_chuncks2(dir1, audio_fname_prefix, fc, vfc );

    marker = 2*ones(size(video_pulse_start_ts,1),1);
    comb_marker_pulses = [video_pulse_start_ts marker];
    
%     pul = figure;
%     plot(tmp,'k')
%     hold on
%     plot(comb_marker_pulses(:,1),comb_marker_pulses(:,2),'r.')

    cd (dir2)
    fname = sprintf('%s.seq',video_fname_prefix);
    fname_path = [dir2 fname];
    strctMovInfo = fnReadSeqInfo_jpn(fname_path);
    clear ts df_ts
    ts = strctMovInfo.m_afTimestamp;
    df_ts = diff(ts);
    dist = figure;
    hist(df_ts,0:0.001:(max(df_ts)+0.01))
%     xlim([(min(df_ts)-0.1) (max(df_ts)+0.1)])
    foo = sprintf('%s_video_pulse_start_ts',video_fname_prefix);
    save(foo,'video_pulse_start_ts')
%     saveas(pul,sprintf('%s_video_pulses.jpg',video_fname_prefix));
    saveas(dist,sprintf('%s_timestamp_distribution.jpg',video_fname_prefix));
else
    cd (dir2)
    foo = sprintf('%s_video_pulse_start_ts',video_fname_prefix);
    load (foo)
    dist = figure('Visible','off');
%     pul = figure('Visible','off');
end

end

