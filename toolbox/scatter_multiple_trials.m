function scatter_multiple_trials(base_dir_name, ...
                                 data_analysis_dir_name, ...
                                 date_str, ...
                                 letter_str, ...
                                 n_vocs_per_trial_max, ...
                                 n_vocs_per_job_max, ...
                                 use_cluster, ...
                                 blob_from_voc_indicators_and_ancillary, ...
                                 args_template_base, ...
                                 verbosity)  %#ok

% base_dir_name a string
% date_str, letter_str each a cell array of strings
m_file_path=mfilename('fullpath');
toolbox_dir_path=fileparts(m_file_path);
project_dir_path=fileparts(toolbox_dir_path);
exe_path=fullfile(project_dir_path,'bin','run_feval_blob_from_voc_indicators_and_ancillary_executable_custom.sh');
mcr_path='/home/taylora/software/MATLAB/R2012a';
%temp_dir_name=fullfile(project_dir_path,'temp');

n_trials=length(date_str);
%blob_per_voc_per_trial=cell(n_trials,1);
for i_trial=1:n_trials
  i_trial  %#ok
  date_str_this_trial=date_str{i_trial};
  letter_str_this_trial=letter_str{i_trial};
  n_vocs_so_far_this_trial=0;
  
  job_inputs_dir_name =fullfile(project_dir_path,sprintf('job_inputs_%s_%s' ,date_str_this_trial,letter_str_this_trial));
  job_outputs_dir_name=fullfile(project_dir_path,sprintf('job_outputs_%s_%s',date_str_this_trial,letter_str_this_trial));
  job_stdouts_dir_name=fullfile(project_dir_path,sprintf('job_stdouts_%s_%s',date_str_this_trial,letter_str_this_trial));
  job_stderrs_dir_name=fullfile(project_dir_path,sprintf('job_stderrs_%s_%s',date_str_this_trial,letter_str_this_trial));

  % load the per-trial ancillary data
  [syl_name,i_start,i_end,f_lo,f_hi, ...
   r_head,r_tail,R,Temp, ...
   dx,x_grid,y_grid,in_cage]= ...
    ssl_trial_overhead(base_dir_name, ...
                       data_analysis_dir_name, ...
                       date_str_this_trial, ...
                       letter_str_this_trial);  %#ok
  
%   % temporary cheat                   
%   i_goal=35145631; margin=10000; keep=(i_goal-margin<=i_start);
%   syl_name=syl_name(keep);
%   i_start=i_start(keep);
%   i_end=i_end(keep);
%   f_lo=f_lo(keep);
%   f_hi=f_hi(keep);
%   r_head=r_head(:,keep,:);
%   r_tail=r_tail(:,keep,:);
                     
%   % store the per-trial data for return
%   per_trial_ancillary(i_trial).date_str=date_str{i_trial};
%   per_trial_ancillary(i_trial).letter_str=letter_str{i_trial};
%   per_trial_ancillary(i_trial).R=R;
%   per_trial_ancillary(i_trial).Temp=Temp;
%   per_trial_ancillary(i_trial).dx=dx;
%   if args_template_base.return_big_things
%     per_trial_ancillary(i_trial).x_grid=x_grid;
%     per_trial_ancillary(i_trial).y_grid=y_grid;
%     per_trial_ancillary(i_trial).in_cage=in_cage;
%   end
  
  % iterate over the vocalizations in this trial
  n_vocs_this_trial=length(syl_name)  %#ok
  n_vocs_to_process_this_trial=min(n_vocs_this_trial,n_vocs_per_trial_max);
  n_jobs_this_trial=ceil(n_vocs_to_process_this_trial/n_vocs_per_job_max);
  n_mice=size(r_head,3);
  args_common=args_template_base;
  args_common.R=R;
  args_common.Temp=Temp;
  args_common.dx=dx;
  %args_template.x_grid=x_grid;
  %args_template.y_grid=y_grid;
  %args_template.in_cage=in_cage;
  args_common.x_grid=[];  % Save space in input file, will re-generate on cluster
  args_common.y_grid=[];
  args_common.in_cage=[];
  tic
  for i_job_this_trial=1:n_jobs_this_trial
    i_voc_this_trial_first=(i_job_this_trial-1)*n_vocs_per_job_max+1;
    if (i_job_this_trial<n_jobs_this_trial)
      i_voc_this_trial_last=(i_job_this_trial-1)*n_vocs_per_job_max+n_vocs_per_job_max;
    else
      i_voc_this_trial_last=n_vocs_to_process_this_trial;
    end
    n_vocs_this_job=i_voc_this_trial_last-i_voc_this_trial_first+1;
    i_voc_this_trial_first  %#ok
    %args=repmat(args_template,[n_vocs_this_job 1]);
    %argses=struct([]);
    syl_name_this=syl_name(i_voc_this_trial_first:i_voc_this_trial_last);
    argses=struct('syl_name',syl_name_this);
    [argses(:).i_start]=csl_from_vector(i_start(i_voc_this_trial_first:i_voc_this_trial_last));
    [argses(:).i_end]=csl_from_vector(i_end(i_voc_this_trial_first:i_voc_this_trial_last));  
    [argses(:).f_lo]=csl_from_vector(f_lo(i_voc_this_trial_first:i_voc_this_trial_last));  
    [argses(:).f_hi]=csl_from_vector(f_hi(i_voc_this_trial_first:i_voc_this_trial_last));
    for j=1:n_vocs_this_job
      i_voc=i_voc_this_trial_first+j-1;
      argses(j).r_head=reshape(r_head(:,i_voc,:),[2 n_mice]);  
      argses(j).r_tail=reshape(r_tail(:,i_voc,:),[2 n_mice]);
    end
    
    % Make sure all the dirs we need exist
    if ~exist(job_inputs_dir_name,'dir')
      mkdir(job_inputs_dir_name);
    end
    if ~exist(job_outputs_dir_name,'dir')
      mkdir(job_outputs_dir_name);
    end
    if ~exist(job_stdouts_dir_name,'dir')
      mkdir(job_stdouts_dir_name);
    end
    if ~exist(job_stderrs_dir_name,'dir')
      mkdir(job_stderrs_dir_name);
    end
    
    input_file_name= ...
      fullfile(job_inputs_dir_name, ...
               sprintf('input_%s_%s_%s.mat',date_str_this_trial,letter_str_this_trial,syl_name_this{1}));
    output_file_name= ...
      fullfile(job_outputs_dir_name, ...
               sprintf('output_%s_%s_%s.mat',date_str_this_trial,letter_str_this_trial,syl_name_this{1}));
    if ~exist(output_file_name,'file') ,
      save(input_file_name,'blob_from_voc_indicators_and_ancillary', ...
                           'base_dir_name', ...
                           'date_str_this_trial', ...
                           'letter_str_this_trial', ...
                           'args_common', ...
                           'argses', ...
                           'verbosity');
      if use_cluster ,                   
        qsub_str=sprintf('qsub -l short=true -A egnorr -b yes -e "%s" -o "%s" "%s" "%s" "%s" "%s"', ...
                         job_stderrs_dir_name, ...
                         job_stdouts_dir_name, ...
                         exe_path, ...
                         mcr_path, ...
                         output_file_name, ...
                         input_file_name);
        fprintf('%s\n',qsub_str);                
        system(qsub_str);         
        pause(0.005);
      else
        feval_blob_from_voc_indicators_and_ancillary_executable(output_file_name,input_file_name);
      end
    end
    
    n_vocs_so_far_this_trial=n_vocs_so_far_this_trial+1;
    if n_vocs_so_far_this_trial>=n_vocs_per_trial_max
      break
    end
  end
  toc
  % blob_per_voc_per_trial{i_trial}=blob_per_voc;
end

end
