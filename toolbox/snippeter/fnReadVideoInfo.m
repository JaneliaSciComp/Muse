function strctVideoInfo = fnReadVideoInfo(strFileName)
% Wrapper to read frames from video files.
% 
%
%Copyright (c) 2008 Shay Ohayon, California Institute of Technology. 
% This file is a part of a free software. you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation (see GPL.txt)

%global g_strVideoWrapper

[dummy,dummy,strExt] = fileparts(strFileName);  %#ok
if strcmpi(strExt,'.seq')
  strctVideoInfo = fnReadSeqInfo(strFileName);
  return;
end;

if strcmpi(strExt,'.avi') || strcmpi(strExt,'.wmv')
  vidInfo = mmfileinfo(strFileName);
  strctVideoInfo.m_strFileName = strFileName;
  strctVideoInfo.m_fFPS = 30;
  strctVideoInfo.m_iNumFrames = floor(vidInfo.Duration*strctVideoInfo.m_fFPS)-1;
  strctVideoInfo.m_iHeight = vidInfo.Video.Height;
  strctVideoInfo.m_iWidth = vidInfo.Video.Width;
  return;
end;

error('fnReadVideoInfo:unsupportedFileType', ...
      'Unable to get info on %s files.', ...
      strExt);