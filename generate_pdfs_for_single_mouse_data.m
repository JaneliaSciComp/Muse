% calculate stuff for all vocalizations
verbosity=0;
args.read_from_map_cache=true;
args.write_to_map_cache=false;
args.pdf_page_dir_name='single_mouse_data_single_voc_pdfs';
args.quantify_confidence=true;  % calculate P-vals, CRs
args.return_big_things=true;  % return the full map and other large
                              % data structures from 
                              % pdf_from_voc_indicators_and_ancillary()

% identifying info for each trial
date_str=cell(0,1);
letter_str=cell(0,1);
date_str{end+1}='06052012';
letter_str{end+1}='D';
date_str{end+1}='06062012';
letter_str{end+1}='E';
date_str{end+1}='06102012';
letter_str{end+1}='E';
date_str{end+1}='06112012';
letter_str{end+1}='D';
date_str{end+1}='06122012';
letter_str{end+1}='D';
date_str{end+1}='06122012';
letter_str{end+1}='E';
date_str{end+1}='06132012';  % this is the one with the least vocs
letter_str{end+1}='D';
date_str{end+1}='06132012';
letter_str{end+1}='E';
n_vocs_per_trial_max=inf;  % max number of vocs to do per trial

base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_sys_test';
data_analysis_dir_name='Data_analysis';

% delete the existing dir, if present
if ~isempty(dir(args.pdf_page_dir_name))
  system(sprintf('rm -rf "%s"',args.pdf_page_dir_name));
end

% make a pdf page for each voc
tic
blob= ...
  map_multiple_trials(base_dir_name, ...
                      data_analysis_dir_name, ...
                      date_str, ...
                      letter_str, ...
                      n_vocs_per_trial_max, ...
                      @pdf_from_voc_indicators_and_ancillary, ...
                      args, ...
                      verbosity);
toc


