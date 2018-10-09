function varargout = tif2mj(input_base,output_base,compression_lvl,delete_file,maxThreads)
%tif2mj Converts raw tiff files into compressed mj2 and DELETES!@! original file.
%
% Inputs:
%   input_base: raw tif folder
%   output_base: output folder of compressed mj2 file
%   compression_lvl: compression level, e.g. 5
%   delete_file: flag to delete (1) the raw file or not (0)
%
% Examples:
%   tif2mj('./data/', './test/','5','0')
%

% $Author: base $	$Date: 2018/10/10 09:28:24 $	$Revision: 0.1 $
% Copyright: HHMI 2016

%% mcc -m -R -nojvm -v <function.m> -d <outfolder/>  -a <addfolder>
if 0
    % numcores = 4;
    compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/tif2mj/tif2mj';
    mkdir(fileparts(compiledfunc))
    unix(sprintf('umask g+rxw %s',fileparts(compiledfunc)))
    sprintf('mcc -m -v -R -singleCompThread %s/tif2mj.m -d %s',pwd,fileparts(compiledfunc))
    % unix(sprintf('mcc -m -v -R -singleCompThread %s/tif2mj.m -d %s',pwd,fileparts(compiledfunc)))
    unix(sprintf('mcc -m -v %s/tif2mj.m -d %s',pwd,fileparts(compiledfunc)))
    unix(sprintf('chmod g+rwx %s',compiledfunc))
end
%%

if nargin<2
    myshfile = '2017-09-25.sh'
    deployment(myshfile)
    return
end
unix('umask 0022')

if nargin<2
    output_base = input_base;
    compression_lvl = '10';
    delete_file = '0';
    maxThreads = '1';
elseif nargin<3
    compression_lvl = '10';
    delete_file = '0';
    maxThreads = '1';
elseif nargin<4
    delete_file = '0';
    maxThreads = '1';
elseif nargin<5
    maxThreads = '1';
end


varargout{1} = 0;
comp = str2num(compression_lvl);
delete_file = str2num(delete_file);
maxThreads = str2num(maxThreads);
LASTN = maxNumCompThreads(maxThreads); % sets the number of threads that will be used

if ~exist(output_base,'dir');mkdir(output_base);end
% check existing tif files
my_tif_files = dir(fullfile(input_base,'*.tif'));
for i_tif = 1:length(my_tif_files)
    input_name = my_tif_files(i_tif).name;
    input_tif_file = fullfile(input_base,input_name);
    [~,outname] = fileparts(input_name);
    outname = sprintf('%s_comp-%02d.mj2',outname,comp);
    output_mj2_file = fullfile(output_base,outname);
    
    % reads tif file
    imgdata = deployedtiffread(input_tif_file);
    v = VideoWriter(output_mj2_file,'Motion JPEG 2000');
    v.CompressionRatio = comp;
    v.MJ2BitDepth = 16;
    open(v)
    writeVideo(v,reshape(imgdata,[size(imgdata,1),size(imgdata,2),1,size(imgdata,3)]))
    close(v)
    unix(sprintf('chmod g+rwx %s',output_mj2_file))

    % delete raw file, -f forces to delete write protected files
    if delete_file; unix(sprintf('rm -f %s',input_tif_file)); end
end
%%
end
%%
function [Iout] = deployedtiffread(fileName,slices)
%DEPLOYEDTIFFREAD Summary of this function goes here
%
% [OUTPUTARGS] = DEPLOYEDTIFFREAD(INPUTARGS) Explain usage here
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/08/21 12:26:16 $	$Revision: 0.1 $
% Copyright: HHMI 2015
warning off
info = imfinfo(fileName, 'tif');
if nargin<2
    slices = 1:length(info);
end
wIm=info(1).Width;
hIm=info(1).Height;
numIm = numel(slices);
Iout  = zeros(hIm, wIm, numIm,'uint16');

for i=1:numIm
    Iout(:,:,i) = imread(fileName,'Index',slices(i),'Info',info);
end

end

%%
function deployment(myshfile,myfilefolder)
%% (comp, tiffile, videofolder, outfolder,outname)
%%
% mcc -m -R -nojvm -v compressbrain_cluster.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_conversion2/ -o compressbrain_cluster
%%
% clear all
% read from a folder
% myshfile
if nargin<2
    myfilefolder = '/groups/mousebrainmicro/mousebrainmicro/data/';
end
opt.myshfile = myshfile;

% myfilefolder = '/groups/mousebrainmicro/mousebrainmicro/from_tier2/data/'
% opt.myshfile = '2015-02-27.sh'
% opt.myshfile = '2015-09-01.sh'
% opt.myshfile = '2014-11-24.sh'
% opt.myshfile = '2015-07-11.sh'
% opt.myshfile = '2016-07-18.sh'
% opt.myshfile = '2016-04-04_rem.sh'
% opt.myshfile = '2016-09-25_r1-1000.sh'
% opt.myshfile = '2016-10-25.sh'
% opt.myshfile = '2016-10-31.sh'
% opt.myshfile = '2016-10-31.sh'
% opt.myshfile = '2016-12-05.sh'
% opt.myshfile = '2017-02-13.sh'
% opt.myshfile = '2017-01-15.sh'
% opt.myshfile = '2017-02-22.sh'
% opt.myshfile = '2017-04-19.sh'
% opt.myshfile = '2016-03-21.sh'
% opt.myshfile = '2016-03-21.sh'
% opt.myshfile = '2016-09-12.sh'
% opt.myshfile = '2014-06-24.sh'

brain = opt.myshfile(1:10);
% check Tiling folder exists
opt.inputfolder = fullfile(myfilefolder,sprintf('%s/Tiling',brain))
if exist(opt.inputfolder,'dir')
else
    opt.inputfolder = fullfile(myfilefolder,sprintf('%s/',brain));
end
    
opt.seqtemp = fullfile(opt.inputfolder,sprintf('filelist.txt'))
unix(sprintf('rm %s',opt.seqtemp))
clear args
args.level = 3;
args.ext = 'tif';
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
if 0 % use array tasks
    %%
    % create main entry
    mkdir(fullfile(opt.inputfolder,'comptasklist'))
    fid = fopen(fullfile(opt.inputfolder,'comp_task_list.sh'),'w');
    fprintf(fid,'#!/bin/sh\n');
    totnumtask=1000;
    range=[1 totnumtask]
    fprintf(fid,'qsub -N comp -t %d-%d -j y -o /dev/null -b y -cwd -V ''%s/comptasklist/comp.${SGE_TASK_ID}.sh''',...
        range(1),range(2),opt.inputfolder);
    fclose(fid);
    unix(sprintf('chmod g+x %s',fullfile(opt.inputfolder,'comp_task_list.sh')));
    
    %%
    comp = 10%[1 5 10 20 40 80];
    delete_file = 1;
    numcores = 3;
    howlong = 3*60;% in seconds
    s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
    %find number of random characters to choose from
    numRands = length(s);
    %specify length of random string to generate
    sLength = 10;

    start = round(linspace(1,numfiles+1,totnumtask+1));
    % create tasks
    for idxtask=1:totnumtask
        % create task
        myfile = sprintf('%s/comptasklist/comp.%d.sh',opt.inputfolder,idxtask);
        fid=fopen(myfile,'w');
        for idx=start(idxtask):start(idxtask+1)-1
            tiffile = myfiles{idx};
            %outfile = tiffile(1:end-4);
            [videofolder,outname] = fileparts(tiffile);
            randString = s( ceil(rand(1,sLength)*numRands) );
            name = sprintf('c-%05d-%s',idx,randString);
            args = sprintf('''/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_conversion2/compressbrain_cluster %d %s %s %s %d> output.log''',comp,tiffile,videofolder,outname,delete_file);
            mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,howlong,name,args);
            fwrite(fid,mysub);
        end
        fclose(fid);
        unix(sprintf('chmod g+x %s',myfile));
    end
    
    %%
else
    %%
    myshfile = opt.myshfile;
    % videofolder = '/tier2/mousebrainmicro/mousebrainmicro/cluster/compressionExperiment/out'
    comp = 10%[1 5 10 20 40 80];
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
        tiffile = myfiles{idx};
        %outfile = tiffile(1:end-4);
        [videofolder,outname] = fileparts(tiffile);
        randString = s( ceil(rand(1,sLength)*numRands) );
        name = sprintf('c-%05d-%s',idx,randString);
        args = sprintf('''/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_conversion2/compressbrain_cluster %d %s %s %s %d''',comp,tiffile,videofolder,outname,delete_file);
%         mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,howlong,name,args);
        mysub = sprintf('bsub -n%d -R"affinity[core(1)]" -We %d -J %s -o /dev/null %s\n',numcores,howlong/60,name,args);
        fwrite(fid,mysub);
    end
    unix(sprintf('chmod +x %s',myshfile));
end
end
