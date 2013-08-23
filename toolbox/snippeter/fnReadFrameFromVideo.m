function Output = fnReadFrameFromVideo(strctVideoInfo, iFrame)
% Wrapper to read frames from video files.
%
%
%Copyright (c) 2008 Shay Ohayon, California Institute of Technology.
% This file is a part of a free software. you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation (see GPL.txt)

%global g_strVideoWrapper
strFileName = strctVideoInfo.m_strFileName;

bForceBW = 1;
[dummy,dummy, strExt] = fileparts(strFileName);  %#ok
if strcmpi(strExt,'.seq')
    Output = fnReadFrameFromSeq(strctVideoInfo, iFrame);
    if bForceBW && size(Output,3) > 1
        Output = Output(:,:,1);
    end;
elseif strcmpi(strExt,'.avi') || strcmpi(strExt,'.wmv')
   vidObj = mmreader(strFileName);
   a3fFrameThis=vidObj.read(iFrame);  % uint8, w x h x 3
   Output=uint8(mean(double(a3fFrameThis),3));  
else
    error('fnReadFrameFromVideo:unsupportedFileType', ...
          'Unable to read %s files.', ...
          strExt);
end

