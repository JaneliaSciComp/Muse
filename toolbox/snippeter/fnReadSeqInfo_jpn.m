function strctMovInfo = fnReadSeqInfo_jpn(strSeqFileName)
% Credits go for Poitr Dollar for the initial code of reading the SEQ
% files.


hFileID = fopen(strSeqFileName);
fseek(hFileID,0,'bof');
% first 4 bytes store OxFEED, next 24 store 'Norpix seq  '
if ~(strcmp(sprintf('%X',fread(hFileID,1,'uint32')),'FEED'))
    % Attempt to fix SEQ header.
    fclose(hFileID);
    fprintf('Header is corrupted for file %s!\n', strSeqFileName);
    strResponse = input('Do you want to fix the file [Y]/[N]? : ','s');
    bFix = strResponse(1) == 'Y' || strResponse(1) == 'y';
    if bFix
        fnFixSeqHeader(strSeqFileName); 
        hFileID = fopen(strSeqFileName);
        fseek(hFileID,0,'bof');
        assert(strcmp(sprintf('%X',fread(hFileID,1,'uint32')),'FEED'));
    else
        strctMovInfo = [];
        return;
    end;
end;
assert(strcmp(char(fread(hFileID,10,'uint16'))','Norpix seq')); %#ok<FREAD>
fseek(hFileID,4,'cof');
% next 8 bytes for version and header size (1024), then 512 for descr
iVersion=fread(hFileID,1,'int32'); 
assert(fread(hFileID,1,'uint32')==1024);
fseek(hFileID,512,'cof');
% read in more strctMovInfo
afBuffer=fread(hFileID,9,'uint32'); 
assert(afBuffer(8)==0);
fps = fread(hFileID,1,'float64');
% store strctMovInformation in strctMovInfo struct
strctMovInfo=struct( 'm_strFileName', strSeqFileName,...
                     'm_iWidth',afBuffer(1), ...
                     'm_iHeight',afBuffer(2), ...
                     'm_iImageBitDepth',afBuffer(3), ...
                     'm_iImageBitDepthReal',afBuffer(4), ...
                     'm_iImageSizeBytes',afBuffer(5), ...
                     'm_iImageFormat',afBuffer(6), ...
                     'm_iNumFrames',afBuffer(7), ...
                     'm_iTrueImageSize', afBuffer(9),...
                     'm_fFPS',fps, ...
                     'm_iSeqiVersion',iVersion);
fclose(hFileID);

% Automatically generate seeking strctMovInfo if not exist
[strPath, strFileName, strExt] = fileparts(strSeqFileName);
if isunix || ismac
    strSeekFilename = [strPath,'/',strFileName,'.mat'];
else
    if length(strPath) == 3 && strPath(3) == '\'
        strSeekFilename = [strPath,strFileName,'.mat'];
    else
        strSeekFilename = [strPath,'\',strFileName,'.mat'];
    end
end;
if ~exist(strSeekFilename,'file')
    [aiSeekPos, afTimestamp] = ...
        fnGenerateSeqSeekInfo(strctMovInfo,strctMovInfo.m_iNumFrames);
    strctMovInfo.m_aiSeekPos = aiSeekPos;
    strctMovInfo.m_afTimestamp = afTimestamp;
    save(strSeekFilename,'strSeqFileName','aiSeekPos','afTimestamp');
else
    strctTmp = load(strSeekFilename);
    strctMovInfo.m_aiSeekPos = strctTmp.aiSeekPos;
    strctMovInfo.m_afTimestamp = strctTmp.afTimestamp;
end;

return;


function [aiSeekPos, afTimestamp] = fnGenerateSeqSeekInfo(strctMovInfo, iNumFrames)
aiSeekPos = zeros(1, iNumFrames);
afTimestamp = zeros(1, iNumFrames);

hFileID = fopen(strctMovInfo.m_strFileName);
fseek(hFileID, aiSeekPos(1),'bof');
if iNumFrames > 1000
       fprintf('Generating seek info for %d frames, please wait...\n', iNumFrames);
end;

if strctMovInfo.m_iImageFormat == 100 || strctMovInfo.m_iImageFormat == 200 
    aiSeekPos = 1024 + (0:iNumFrames)*strctMovInfo.m_iTrueImageSize;
    % Read timestamp info...
    for iIter = 0:iNumFrames-1
        fseek(hFileID,1024+iIter*strctMovInfo.m_iTrueImageSize+strctMovInfo.m_iImageSizeBytes,'bof');
        iA = fread(hFileID,1,'uint32');
        iB = fread(hFileID,1,'uint16');
        afTimestamp(iIter+1) =  double(iA)+ double(iB)/1000;
    end
elseif  strctMovInfo.m_iImageFormat == 102 || strctMovInfo.m_iImageFormat == 201
    aiSeekPos(1) = 1024;
    if fnNewSEQFileType(strctMovInfo.m_strFileName)
        fprintf('New SEQ file type detected!\n');
        iOffset = 8;
    else
        iOffset = 16;
    end
    % Compressed seq
    for iIter = 0:iNumFrames-1
        if mod(iIter,10000) == 0
            fprintf('Passed frame %d\n',iIter);
        end;
        fseek(hFileID,aiSeekPos(iIter+1),'bof');
        iCurrImageSizeBytes = fread(hFileID,1,'uint32');
        fseek(hFileID,iCurrImageSizeBytes-4,'cof');
        iA = fread(hFileID,1,'uint32');
        iB = fread(hFileID,1,'uint16');
        afTimestamp(iIter+1) =  double(iA)+ double(iB)/1000;
        if iIter ~= iNumFrames-1
            aiSeekPos(iIter+2) = aiSeekPos(iIter+1) + iCurrImageSizeBytes + iOffset;
        end
    end
        
    aiSeekPos = aiSeekPos + 4; % Skip image size and go directly to the JPG header.
else
    assert(false);
end;
fprintf('Done!\n');
fclose(hFileID);
return;


function bNewFileType = fnNewSEQFileType(strFileName)
% The direct method.. try to parse the JPG header.
% if it fails, it means we jumped too much!
hFileID = fopen(strFileName);
fseek(hFileID,1024,'bof');
iFirstFrameSizeInBytes = fread(hFileID,1,'uint32');
fseek(hFileID,iFirstFrameSizeInBytes-4,'cof');
iADummy = fread(hFileID,1,'uint32'); % Read time stamp...
iBDummy = fread(hFileID,1,'uint16'); % Read time stamp...
iNextPositionNewFile = 1028 +iFirstFrameSizeInBytes + 8;
try
    X = parsejpg8(strFileName, iNextPositionNewFile);
    bNewFileType = true;
catch
    bNewFileType = false;
end

fclose(hFileID);
return;
