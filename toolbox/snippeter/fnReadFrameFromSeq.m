function a2iFrame = fnReadFrameFromSeq(strctMovInfo, iFrame)

if strctMovInfo.m_iImageFormat == 100 || strctMovInfo.m_iImageFormat == 200
  % Read uncompressed frame from video
  hFileID = fopen(strctMovInfo.m_strFileName);
  fseek(hFileID, strctMovInfo.m_aiSeekPos(iFrame), 'bof');
  dataRaw=fread(hFileID,strctMovInfo.m_iImageSizeBytes,'uint8=>uint8');
  a2iFrame = reshape(dataRaw, ...
                     strctMovInfo.m_iWidth, ...
                     strctMovInfo.m_iHeight)';
  fclose(hFileID);
elseif strctMovInfo.m_iImageFormat == 102 || strctMovInfo.m_iImageFormat == 201
  a2iFrame = parsejpg8(strctMovInfo.m_strFileName, ...
                       strctMovInfo.m_aiSeekPos(iFrame));
  if size(a2iFrame,3) == 3
    a2iFrame = rgb2gray(a2iFrame);
  end  
else
  assert(false);
end

end
