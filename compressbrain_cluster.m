function compressbrain_cluster(comp, tiffile, videofolder,outname,delete_file)
%compressbrain Converts raw tiff files into compressed mj2 and DELETES!@! original file.
% This is based on runvonv.m
%
% RUNCONV(comp, tiffile, videofolder, outfolder,outname)
%
% Inputs:
%   comp: compression level, e.g. 5
%   tiffile: raw tif location
%   videofolder: location of the output mj2 file
%   outname: prefix for output file
%   delete_file: flag to delete (1) the raw file or not (0)
%
% Examples:
%   compressbrain_cluster('5', './data/I2.tif', './test','myvideo','0')
%
% See also: runconv.m

% $Author: base $	$Date: 2016/01/21 09:28:24 $	$Revision: 0.1 $
% Copyright: HHMI 2016
%%
unix('umask 0022')
comp = str2num(comp);
delete_file = str2num(delete_file);
mkdir(videofolder)
videoname = fullfile(videofolder,sprintf('%s_comp-%02d.mj2',outname,comp));
%%
% reads tif file
imgdata = deployedtiffread(tiffile);
v = VideoWriter(videoname,'Motion JPEG 2000');
v.CompressionRatio = comp;
v.MJ2BitDepth = 16;
open(v)
writeVideo(v,reshape(imgdata,[size(imgdata,1),size(imgdata,2),1,size(imgdata,3)]))
close(v)

% delete raw file, -f forces to delete write protected files
%     unix(sprintf('rm -f %s',tiffile))
if delete_file
    unix(sprintf('rm -f %s',tiffile))
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
function deployment
%% (comp, tiffile, videofolder, outfolder,outname)
%%
% mcc -m -R -nojvm -v compressbrain_cluster.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_conversion/ -o compressbrain_cluster
%%
clear all
% read from a folder
opt.myshfile = '2015-02-27.sh'
opt.myshfile = '2015-09-01.sh'
opt.myshfile = '2014-11-24.sh'
opt.myshfile = '2015-07-11.sh'
opt.myshfile = '2016-07-18.sh'
opt.myshfile = '2016-04-04_rem.sh'
opt.myshfile = '2016-09-25_r1-1000.sh'
brain = opt.myshfile(1:10);
opt.inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/from_tier2/data/%s/Tiling',brain)
opt.seqtemp = fullfile(opt.inputfolder,sprintf('filelist_range1-1413.txt'))
unix(sprintf('rm %s',opt.seqtemp))
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
            args = sprintf('''/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_conversion/compressbrain_cluster %d %s %s %s %d> output.log''',comp,tiffile,videofolder,outname,delete_file);
            mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,howlong,name,args);
            fwrite(fid,mysub);
        end
        fclose(fid);
        unix(sprintf('chmod g+x %s',myfile));
        
    end
    %%
else
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
        args = sprintf('''/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/compiled/compiledfiles_conversion/compressbrain_cluster %d %s %s %s %d> output.log''',comp,tiffile,videofolder,outname,delete_file);
        mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o ~/logs -b y -cwd -V %s\n',numcores,howlong,name,args);
        fwrite(fid,mysub);
    end
    unix(sprintf('chmod +x %s',myshfile));
end
end


% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %% read from yml
% clear all
% close all
% clc
% if 1
% %     load frames-2015-06-19-johan-full
%     load ./data/frames-2015-04-24_orig frames
% else
%     frames = LoadBasicYAML('/tier2/mousebrainmicro/mousebrainmicro/data/2015-04-24/Tiling/tilebase.cache.orig.yml');
%     frames = UnionBBoxes(frames);
%     save ./data/frames-2015-04-24_orig frames
% end

% myfile = sprintf('./shfiles/cluster_datazip.sh');
% idxTiles = 1:frames.N;%[2648    3392    3752    5666    7362    8368    8369    8716    8717];

% %%
% % videofolder = '/tier2/mousebrainmicro/mousebrainmicro/cluster/compressionExperiment/out'
% delete_file = 1;
% numcores = 3;
% s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
% %find number of random characters to choose from
% numRands = length(s);
% %specify length of random string to generate
% sLength = 10;
% %-o /dev/null
% fid = fopen(myfile,'w');
% for idx = 1:length(idxTiles)
%     datafolder = fullfile(frames.root,frames.path{idxTiles(idx)});
%     videofolder = datafolder; % write to input folder
%     myfiles = dir([datafolder,'/*.tif']);
%     for fi = 1:2
%         tiffile = fullfile(datafolder, myfiles(fi).name);
%         for comp = 10%[1 5 10 20 40 80];
%             outname = myfiles(fi).name(1:end-4); % use inputname as output
%             % sprintf('IM-%04d_idx-%05d_ch-%d',idx,idxTiles(idx),fi-1);

%             randString = s( ceil(rand(1,sLength)*numRands) );
%             name = sprintf('skel_%05d-%s',idx,randString);
%             args = sprintf('''./compiledfiles_conversion/compressbrain_cluster %d %s %s %s %d> output.log''',comp,tiffile,videofolder,outname,delete_file);
%             mysub = sprintf('qsub -pe batch %d -l short=true -N %s -j y -o ~/logs -b y -cwd -V %s\n',numcores,name,args);
%             fwrite(fid,mysub);
%         end
%     end
% end
% unix(sprintf('chmod +x %s',myfile));
% %%
% end
% function deployment_ilastik
% %% (comp, tiffile, videofolder, outfolder,outname)
% clear all
% close all
% clc

% functionname = ''
% numcores = 8;
% memsize = numcores*7.5*1000;

% load frames-2015-06-19-johan-full
% idxTiles = [2648 3392 3752 5666 7362 8368 8369 8716 8717];

% infolder = '/tier2/mousebrainmicro/mousebrainmicro/cluster/compressionExperiment/out/h5/'
% logfolder = '/nobackup/mousebrainmicro/cluster/2015-06-19-GN3-Takako/logs_comp/'
% mkdir(logfolder)
% out = '/tier2/mousebrainmicro/mousebrainmicro/cluster/compressionExperiment/probs/'
% mkdir(out)

% myfile = sprintf('cluster_conversion_exp_ilp.sh');
% s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
% %find number of random characters to choose from
% numRands = length(s);
% %specify length of random string to generate
% sLength = 10;
% %-o /dev/null
% ilastikloc = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/ilastik-1.1.8.post1-Linux/run_ilastik.sh'
% fid = fopen(myfile,'w');
% for idx = 1:length(idxTiles)
%     for fi = 1:2
%         for comp = [1 5 10 20 40 80];
%             if comp~=1
%                 ilpfile = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/compressed_data_for_training/compressionTest_cmp%02d.ilp',comp);
%             else
%                 ilpfile = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/GN1.ilp';
%             end
%             outname = sprintf('IM-%04d_idx-%05d_ch-%d_comp-%02d',idx,idxTiles(idx),fi-1,comp);
%             infiles = fullfile(infolder,outname);

%             randString = s( ceil(rand(1,sLength)*numRands) );
%             name = sprintf('exp_%05d-%s',idx,randString);
%             logfile=fullfile(logfolder,sprintf('ilastikclassifierrun_%05d-%s.txt',idx,randString));
%             outputformat=sprintf('''%s_Probabilities.h5''',fullfile(out,outname));
%             args = sprintf('''%s --headless --logfile=%s --project=%s --output_format=hdf5 --output_filename_format=%s %s''',...
%                 ilastikloc,logfile,ilpfile,outputformat,[infiles,'.h5']);
%             mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d qsub -pe batch %d -l short=true -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,memsize,numcores,name,args);
%             fwrite(fid,mysub);
%         end
%     end
% end
% unix(sprintf('chmod +x %s',myfile));
% end

% function deployment_skel
% %%
% numcores = 8;
% threshold='"[.5,.9]"';
% % threshold='"[.3,.5,.7,.9]"'

% switch 'GN3'
%     case 'GN1'
%         out = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/GN1_comp10_skel'
%         infolder = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/GN1_comp10_probs/'
%         myfile = 'cluster_skel_GN1_comp10.sh';
%     case 'GN3'
%         out = '/tier2/mousebrainmicro/mousebrainmicro/cluster/compressionExperiment/skel'
%         infolder = '/tier2/mousebrainmicro/mousebrainmicro/cluster/compressionExperiment/probs/'
%         myfile = 'cluster_skel_GN3_compexps.sh';
% end

% mkdir(out)
% s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';

% %find number of random characters to choose from
% numRands = length(s);
% %specify length of random string to generate
% sLength = 10;

% myfiles = dir([infolder,'*.h5']);
% %-o /dev/null
% fid = fopen(myfile,'w');
% for ii=1:length(myfiles)
%     %%
%     %generate random string
%     inputimage = fullfile(infolder,myfiles(ii).name);
%     randString = s( ceil(rand(1,sLength)*numRands) );
%     name = sprintf('skel_%05d-%s',ii,randString);
%     args = sprintf('''./compiledfiles_skel/ilp_bpTile %s %s %s> output.log''',inputimage,threshold,out);
%     mysub = sprintf('qsub -pe batch %d -l short=true -N %s -j y -o ~/logs -b y -cwd -V %s\n',numcores,name,args);
%     fwrite(fid,mysub);
% end
% unix(sprintf('chmod +x %s',myfile));



% end