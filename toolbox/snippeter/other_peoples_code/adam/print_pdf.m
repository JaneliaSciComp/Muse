function print_pdf(fig_h,basename)

% prints to a PDF file.  Needs ghostscript to work

tic
if nargin<1 || isempty(fig_h)
  fig_h=gcf;
end
if nargin<2 || isempty(basename)
  if isempty(get(fig_h,'name'))
    basename=sprintf('fig-%03d',fig_h);
  else
    basename=get(fig_h,'name');
  end
end
fprintf(1,'Writing PDF file %s.pdf:\n',basename);
temp_file_path=[tempname '.eps'];
%print(fig_h,'-depsc2','-loose','-adobecset',temp_file_path);
old_vals=set_figure_to_wysiwyg_printing(fig_h);
print(fig_h,'-depsc2','-loose',temp_file_path);
unset_figure_from_wysiwyg_printing(fig_h,old_vals);
if ispc
  if strcmp(computer('arch'),'win64')
    command_name='gswin64c';
  else
    command_name='gswin32c';
  end
else
  command_name='gs';
end
eval_me=...
  sprintf(['! %s -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=pdfwrite ' ...
           '-dAutoRotatePages=/None -sOutputFile="%s.pdf" "%s"'],...
          command_name,...
          basename,...
          temp_file_path);
%eval_me=sprintf('! acrodist /n /q "%s\\%s.eps"',pwd,basename);
eval(eval_me);
delete(temp_file_path);
t=toc;
fprintf(1,'Elapsed time: %0.1f sec\n',t);
