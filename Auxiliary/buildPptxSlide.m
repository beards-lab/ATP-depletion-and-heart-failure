function outFile = buildPptxSlide(outFile, slide)
%BUILDPPTXSLIDE  Write a single-slide .pptx from images and text boxes.
%
%   OUT = BUILDPPTXSLIDE(OUTFILE, SLIDE)
%
%   A .pptx is a ZIP of OOXML parts, so this writes the parts and zips them --
%   no MATLAB Report Generator, no Office automation, no Python. Deliberately
%   minimal: one blank layout, one slide, pictures and text boxes only. That is
%   all a figure slide needs, and it keeps the XML small enough to read.
%
%   SLIDE fields (all optional except .images):
%     .width_in, .height_in   slide size in inches (default 13.3333 x 7.5, 16:9)
%     .bg                     background hex, no '#' (default 'FFFFFF')
%     .images     struct array: .file .x .y .w .h        (inches, from top-left)
%     .texts      struct array: .str .x .y .w .h and optionally
%                              .size (pt, default 18) .bold (logical)
%                              .color (hex, default '1F2933') .align l|ctr|r
%
%   Positions are inches from the top-left corner; OOXML uses EMU internally
%   (914400 per inch) and the conversion happens here.
%
%   Example:
%       s.images(1) = struct('file','a.png','x',0.5,'y',1.2,'w',4.3,'h',2.87);
%       s.texts(1)  = struct('str','Title','x',0.5,'y',0.2,'w',12,'h',0.6, ...
%                            'size',30,'bold',true,'color','1F2933','align','l');
%       buildPptxSlide('out.pptx', s);
%
%   Verified to open in PowerPoint 2021 and LibreOffice Impress.

    EMU = 914400;
    if ~isfield(slide,'width_in')  || isempty(slide.width_in);  slide.width_in  = 13.3333; end
    if ~isfield(slide,'height_in') || isempty(slide.height_in); slide.height_in = 7.5;     end
    if ~isfield(slide,'bg') || isempty(slide.bg); slide.bg = 'FFFFFF'; end
    if ~isfield(slide,'texts'); slide.texts = struct([]); end
    cx = round(slide.width_in*EMU); cy = round(slide.height_in*EMU);

    work = fullfile(tempname); mkdir(work);
    cleanup = onCleanup(@() rmdir(work,'s'));
    mk = @(varargin) iMkdir(work, varargin{:});
    mk('_rels'); mk('ppt','_rels'); mk('ppt','slides','_rels');
    mk('ppt','slideLayouts','_rels'); mk('ppt','slideMasters','_rels');
    mk('ppt','theme'); mk('ppt','media');

    %% ---- media -----------------------------------------------------------
    nI = numel(slide.images);
    for i = 1:nI
        if ~isfile(slide.images(i).file)
            error('buildPptxSlide:noImage','image not found: %s', slide.images(i).file);
        end
        [~,~,e] = fileparts(slide.images(i).file);
        if ~strcmpi(e,'.png')
            error('buildPptxSlide:pngOnly','only .png is supported (got %s)', e);
        end
        copyfile(slide.images(i).file, fullfile(work,'ppt','media',sprintf('image%d.png',i)));
    end

    %% ---- [Content_Types].xml --------------------------------------------
    iWrite(fullfile(work,'[Content_Types].xml'), [ ...
      iHdr() ...
      '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">' ...
      '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' ...
      '<Default Extension="xml" ContentType="application/xml"/>' ...
      '<Default Extension="png" ContentType="image/png"/>' ...
      '<Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/>' ...
      '<Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/>' ...
      '<Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml"/>' ...
      '<Override PartName="/ppt/slides/slide1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>' ...
      '<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>' ...
      '</Types>']);

    %% ---- package rels ----------------------------------------------------
    iWrite(fullfile(work,'_rels','.rels'), [iHdr() iRels({ ...
        'rId1','officeDocument','ppt/presentation.xml'})]);

    iWrite(fullfile(work,'ppt','_rels','presentation.xml.rels'), [iHdr() iRels({ ...
        'rId1','slideMaster','slideMasters/slideMaster1.xml'; ...
        'rId2','slide',      'slides/slide1.xml'; ...
        'rId3','theme',      'theme/theme1.xml'})]);

    iWrite(fullfile(work,'ppt','slideMasters','_rels','slideMaster1.xml.rels'), [iHdr() iRels({ ...
        'rId1','slideLayout','../slideLayouts/slideLayout1.xml'; ...
        'rId2','theme',      '../theme/theme1.xml'})]);

    iWrite(fullfile(work,'ppt','slideLayouts','_rels','slideLayout1.xml.rels'), [iHdr() iRels({ ...
        'rId1','slideMaster','../slideMasters/slideMaster1.xml'})]);

    R = cell(nI,3);
    for i = 1:nI
        R(i,:) = {sprintf('rId%d',i), 'image', sprintf('../media/image%d.png',i)};
    end
    if nI == 0; R = cell(0,3); end
    iWrite(fullfile(work,'ppt','slides','_rels','slide1.xml.rels'), [iHdr() iRels(R)]);

    %% ---- presentation.xml ------------------------------------------------
    iWrite(fullfile(work,'ppt','presentation.xml'), [iHdr() sprintf([ ...
      '<p:presentation xmlns:a="%s" xmlns:r="%s" xmlns:p="%s">' ...
      '<p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId1"/></p:sldMasterIdLst>' ...
      '<p:sldIdLst><p:sldId id="256" r:id="rId2"/></p:sldIdLst>' ...
      '<p:sldSz cx="%d" cy="%d"/><p:notesSz cx="%d" cy="%d"/>' ...
      '</p:presentation>'], iNS('a'), iNS('r'), iNS('p'), cx, cy, cy, cx)]);

    %% ---- theme (minimal but schema-complete) -----------------------------
    iWrite(fullfile(work,'ppt','theme','theme1.xml'), iTheme());

    %% ---- slide master + layout ------------------------------------------
    empt = ['<p:cSld><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/>' ...
            '<p:nvPr/></p:nvGrpSpPr><p:grpSpPr/></p:spTree></p:cSld>'];
    clrMap = ['<p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" ' ...
              'accent2="accent2" accent3="accent3" accent4="accent4" accent5="accent5" ' ...
              'accent6="accent6" hlink="hlink" folHlink="folHlink"/>'];
    iWrite(fullfile(work,'ppt','slideMasters','slideMaster1.xml'), [iHdr() sprintf([ ...
      '<p:sldMaster xmlns:a="%s" xmlns:r="%s" xmlns:p="%s">%s%s' ...
      '<p:sldLayoutIdLst><p:sldLayoutId id="2147483649" r:id="rId1"/></p:sldLayoutIdLst>' ...
      '</p:sldMaster>'], iNS('a'), iNS('r'), iNS('p'), empt, clrMap)]);

    iWrite(fullfile(work,'ppt','slideLayouts','slideLayout1.xml'), [iHdr() sprintf([ ...
      '<p:sldLayout xmlns:a="%s" xmlns:r="%s" xmlns:p="%s" type="blank" preserve="1">%s' ...
      '</p:sldLayout>'], iNS('a'), iNS('r'), iNS('p'), empt)]);

    %% ---- the slide -------------------------------------------------------
    sp = ''; id = 2;
    for i = 1:numel(slide.texts)
        t = slide.texts(i);
        if ~isfield(t,'size')  || isempty(t.size);  t.size  = 18;        end
        if ~isfield(t,'bold')  || isempty(t.bold);  t.bold  = false;     end
        if ~isfield(t,'color') || isempty(t.color); t.color = '1F2933';  end
        if ~isfield(t,'align') || isempty(t.align); t.align = 'l';       end
        sp = [sp sprintf([ ...
          '<p:sp><p:nvSpPr><p:cNvPr id="%d" name="Text %d"/>' ...
          '<p:cNvSpPr txBox="1"/><p:nvPr/></p:nvSpPr>' ...
          '<p:spPr><a:xfrm><a:off x="%d" y="%d"/><a:ext cx="%d" cy="%d"/></a:xfrm>' ...
          '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>' ...
          '<p:txBody><a:bodyPr wrap="square" lIns="0" tIns="0" rIns="0" bIns="0">' ...
          '<a:spAutoFit/></a:bodyPr><a:lstStyle/>%s</p:txBody></p:sp>'], ...
          id, id, round(t.x*EMU), round(t.y*EMU), round(t.w*EMU), round(t.h*EMU), ...
          iParas(t))];                                                  %#ok<AGROW>
        id = id + 1;
    end
    for i = 1:nI
        im = slide.images(i);
        sp = [sp sprintf([ ...
          '<p:pic><p:nvPicPr><p:cNvPr id="%d" name="Picture %d"/>' ...
          '<p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr><p:nvPr/></p:nvPicPr>' ...
          '<p:blipFill><a:blip r:embed="rId%d"/><a:stretch><a:fillRect/></a:stretch></p:blipFill>' ...
          '<p:spPr><a:xfrm><a:off x="%d" y="%d"/><a:ext cx="%d" cy="%d"/></a:xfrm>' ...
          '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>'], ...
          id, id, i, round(im.x*EMU), round(im.y*EMU), round(im.w*EMU), round(im.h*EMU))]; %#ok<AGROW>
        id = id + 1;
    end

    iWrite(fullfile(work,'ppt','slides','slide1.xml'), [iHdr() sprintf([ ...
      '<p:sld xmlns:a="%s" xmlns:r="%s" xmlns:p="%s"><p:cSld>' ...
      '<p:bg><p:bgPr><a:solidFill><a:srgbClr val="%s"/></a:solidFill>' ...
      '<a:effectLst/></p:bgPr></p:bg>' ...
      '<p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>' ...
      '<p:grpSpPr/>%s</p:spTree></p:cSld>%s</p:sld>'], ...
      iNS('a'), iNS('r'), iNS('p'), slide.bg, sp, clrMap)]);

    %% ---- zip -------------------------------------------------------------
    if isfile(outFile); delete(outFile); end
    d = fileparts(outFile); if ~isempty(d) && ~isfolder(d); mkdir(d); end
    parts = [{'[Content_Types].xml','_rels'}, {'ppt'}];
    zip(outFile, parts, work);
    % MATLAB's zip appends .zip when the name has no recognised extension;
    % .pptx is preserved, but check rather than assume.
    if ~isfile(outFile) && isfile([outFile '.zip'])
        movefile([outFile '.zip'], outFile);
    end
    di = dir(outFile);
    fprintf('buildPptxSlide: %s  (%.2f MB, %d images, %d text boxes, %.2f x %.2f in)\n', ...
            outFile, di.bytes/1e6, nI, numel(slide.texts), slide.width_in, slide.height_in);
end

%% ======================= helpers ========================================
function iMkdir(root, varargin); p = fullfile(root, varargin{:});
    if ~isfolder(p); mkdir(p); end
end

function s = iHdr(); s = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'; end

function u = iNS(k)
    switch k
        case 'a'; u = 'http://schemas.openxmlformats.org/drawingml/2006/main';
        case 'r'; u = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships';
        case 'p'; u = 'http://schemas.openxmlformats.org/presentationml/2006/main';
    end
end

function s = iRels(R)
    s = '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">';
    for i = 1:size(R,1)
        s = [s sprintf('<Relationship Id="%s" Type="%s/%s" Target="%s"/>', ...
             R{i,1}, iNS('r'), R{i,2}, R{i,3})];                        %#ok<AGROW>
    end
    s = [s '</Relationships>'];
end

function s = iParas(t)
% One <a:p> per newline-separated line, so a multi-line caption keeps its breaks.
    lines = strsplit(char(t.str), newline);
    s = '';
    for i = 1:numel(lines)
        b = '0'; if t.bold; b = '1'; end
        s = [s sprintf([ ...
            '<a:p><a:pPr algn="%s"/><a:r><a:rPr lang="en-US" sz="%d" b="%s" dirty="0">' ...
            '<a:solidFill><a:srgbClr val="%s"/></a:solidFill>' ...
            '<a:latin typeface="Calibri"/></a:rPr>' ...
            '<a:t>%s</a:t></a:r></a:p>'], ...
            t.align, round(t.size*100), b, t.color, iEsc(lines{i}))];   %#ok<AGROW>
    end
end

function s = iEsc(s)
    s = strrep(s, '&', '&amp;'); s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');  s = strrep(s, '"', '&quot;');
end

function s = iTheme()
% Minimal theme: PowerPoint requires clrScheme, fontScheme and fmtScheme to be
% present and complete, even for a deck that references none of them.
    dk = {'000000','44546A'}; lt = {'FFFFFF','E7E6E6'};
    acc = {'4472C4','ED7D31','A5A5A5','FFC000','5B9BD5','70AD47'};
    c = sprintf('<a:dk1><a:srgbClr val="%s"/></a:dk1><a:lt1><a:srgbClr val="%s"/></a:lt1>', dk{1}, lt{1});
    c = [c sprintf('<a:dk2><a:srgbClr val="%s"/></a:dk2><a:lt2><a:srgbClr val="%s"/></a:lt2>', dk{2}, lt{2})];
    for i = 1:6
        c = [c sprintf('<a:accent%d><a:srgbClr val="%s"/></a:accent%d>', i, acc{i}, i)]; %#ok<AGROW>
    end
    c = [c '<a:hlink><a:srgbClr val="0563C1"/></a:hlink>' ...
           '<a:folHlink><a:srgbClr val="954F72"/></a:folHlink>'];
    fs = ['<a:fontScheme name="Office">' ...
          '<a:majorFont><a:latin typeface="Calibri Light"/><a:ea typeface=""/><a:cs typeface=""/></a:majorFont>' ...
          '<a:minorFont><a:latin typeface="Calibri"/><a:ea typeface=""/><a:cs typeface=""/></a:minorFont>' ...
          '</a:fontScheme>'];
    fill = '<a:solidFill><a:schemeClr val="phClr"/></a:solidFill>';
    ln   = ['<a:ln w="6350" cap="flat" cmpd="sng" algn="ctr">' fill ...
            '<a:prstDash val="solid"/></a:ln>'];
    fmt = ['<a:fmtScheme name="Office">' ...
           '<a:fillStyleLst>' fill fill fill '</a:fillStyleLst>' ...
           '<a:lnStyleLst>' ln ln ln '</a:lnStyleLst>' ...
           '<a:effectStyleLst>' ...
           '<a:effectStyle><a:effectLst/></a:effectStyle>' ...
           '<a:effectStyle><a:effectLst/></a:effectStyle>' ...
           '<a:effectStyle><a:effectLst/></a:effectStyle>' ...
           '</a:effectStyleLst>' ...
           '<a:bgFillStyleLst>' fill fill fill '</a:bgFillStyleLst>' ...
           '</a:fmtScheme>'];
    s = [iHdr() sprintf([ ...
        '<a:theme xmlns:a="%s" name="Office Theme"><a:themeElements>' ...
        '<a:clrScheme name="Office">%s</a:clrScheme>%s%s' ...
        '</a:themeElements><a:objectDefaults/><a:extraClrSchemeLst/></a:theme>'], ...
        iNS('a'), c, fs, fmt)];
end

function iWrite(f, s)
    fid = fopen(f, 'w', 'n', 'UTF-8');
    if fid < 0; error('buildPptxSlide:write','cannot write %s', f); end
    fwrite(fid, unicode2native(s, 'UTF-8')); fclose(fid);
end
