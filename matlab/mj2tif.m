function varargout = mj2tif(input_base,output_base,delete_file,maxThreads)
%tif2mj Converts raw tiff files into compressed mj2 and DELETES!@! original file.
%
% Inputs:
%   input_base: mj2 folder
%   output_base: output folder of decompressed tif file
%   delete_file: flag to delete (1) the raw file or not (0)
%
% Examples:
%   mj2tif('./data/mj2', './data/out','0')
%

% $Author: base $	$Date: 2018/11/20 09:28:24 $	$Revision: 0.1 $
% Copyright: HHMI 2018

%% mcc -m -R -nojvm -v <function.m> -d <outfolder/>  -a <addfolder>
if 0
    % numcores = 4;
    compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/mj2tif/mj2tif';
    mkdir(fileparts(compiledfunc))
    unix(sprintf('umask g+rxw %s',fileparts(compiledfunc)))
    sprintf('mcc -m -v -R -singleCompThread %s/mj2tif.m -d %s',pwd,fileparts(compiledfunc))
    % unix(sprintf('mcc -m -v -R -singleCompThread %s/tif2mj.m -d %s',pwd,fileparts(compiledfunc)))
    unix(sprintf('mcc -m -v %s/mj2tif.m -d %s',pwd,fileparts(compiledfunc)))
    unix(sprintf('chmod g+rwx %s',compiledfunc))
end
%%

if nargin<2
    myshfile = '2017-09-19_decomp.sh'
    deployment(myshfile)
    return
end
unix('umask 0022')

if nargin<2
    output_base = input_base;
    delete_file = '0';
    maxThreads = '1';
elseif nargin<3
    delete_file = '0';
    maxThreads = '1';
elseif nargin<4
    maxThreads = '1';
end


varargout{1} = 0;
delete_file = str2double(delete_file);
maxThreads = str2double(maxThreads);
maxNumCompThreads(maxThreads); % sets the number of threads that will be used

if ~exist(output_base,'dir');mkdir(output_base);end
% check existing tif files
my_mj2_files = dir(fullfile(input_base,'*.mj2'));
for i_mj2 = 1:length(my_mj2_files)
    input_name = my_mj2_files(i_mj2).name;
    input_mj2_file = fullfile(input_base,input_name);
    [~,outname] = fileparts(input_name);
    outname = sprintf('%s_comp-%02d.mj2',outname,comp);
    output_mj2_file = fullfile(output_base,outname);
    
    % reads tif file
    imgdata = deployedtiffread(input_mj2_file);
    v = VideoWriter(output_mj2_file,'Motion JPEG 2000'); %#ok<TNMLP>
    v.CompressionRatio = comp;
    v.MJ2BitDepth = 16;
    open(v)
    writeVideo(v,reshape(imgdata,[size(imgdata,1),size(imgdata,2),1,size(imgdata,3)]))
    close(v)
    try unix(sprintf('chmod g+rwx %s',output_mj2_file));catch;end

    % delete raw file, -f forces to delete write protected files
    if delete_file; unix(sprintf('rm -f %s',input_mj2_file)); end
end
%%
end
%%
function mytiffwrite(outputFileName,imgdata)
imwrite(imgdata(:, :, 1), outputFileName, 'WriteMode', 'overwrite',  'Compression','none');
for K=2:length(imgdata(1, 1, :))
    imwrite(imgdata(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
end
end

%%
function deployment(myshfile,myinputfolder)
%% (comp, tiffile, videofolder, outfolder,outname)
%%
% mcc -m -R -nojvm -v /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/migrated2github/compression/matlab/decompressbrain_cluster.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/decompress

%%
% clear all
% read from a folder
% myshfile
if nargin<2
    myinputfolder = '/nrs/mouselight/Users/base/';
end
opt.myshfile = myshfile;

brain = opt.myshfile(1:10);
% check Tiling folder exists
opt.inputfolder = fullfile(myinputfolder,sprintf('%s/Tiling',brain))
if exist(opt.inputfolder,'dir')
else
    opt.inputfolder = fullfile(myinputfolder,sprintf('%s/',brain));
end
    
opt.seqtemp = fullfile(opt.inputfolder,sprintf('filelist.txt'))
unix(sprintf('rm %s',opt.seqtemp))
clear args
args.level = 3;
args.ext = 'mj2';
if exist(opt.seqtemp, 'file') == 2
    % load file directly
else
    args.fid = fopen(opt.seqtemp,'w');
    recdir(opt.inputfolder,args)
end
fid=fopen(opt.seqtemp,'r');
myfiles = textscan(fid,'%s');
myfiles = myfiles{1};
numfiles = size(myfiles,1);
fclose(fid);

%%
myshfile = opt.myshfile;
delete_file = 1;
numcores = 3;
howlong = 3*60;% in seconds
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
%find number of random characters to choose from
numRands = length(s);
%specify length of random string to generate
sLength = 10;
%-o /dev/null
fid = fopen(myshfile,'w');
for idx = 1:numfiles
    inputfile = myfiles{idx};
    %outfile = tiffile(1:end-4);
    randString = s( ceil(rand(1,sLength)*numRands) );
    name = sprintf('c-%05d-%s',idx,randString);
    args = sprintf('''/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/decompress/decompressbrain_cluster %s %d''',inputfile,delete_file);
    

    % mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,howlong,name,args);
    mysub = sprintf('bsub -n%d -R"affinity[core(1)]" -We %d -J %s -o /dev/null %s\n',numcores,howlong/60,name,args);
    fwrite(fid,mysub);
end
unix(sprintf('chmod +x %s',myshfile));
end
