function decompressbrain_cluster(videofile,delete_file)
%decompressbrain Converts mj2 files into tif and DELETES!@! original file.
%%

if nargin<2
    myshfile = '2017-09-19_decomp.sh'
    deployment(myshfile)
    return
end
unix('umask 0022')

[outfolder,outname] = fileparts(videofile);
parts = strsplit(outname,'.');
outname = [parts{1},'.',parts{2}(1),'.tif'];
outtifname = fullfile(outfolder,outname);

%%
% reads mj2 file
vidObj = VideoReader(videofile);
numFrames = ceil(vidObj.FrameRate*vidObj.Duration);
outimage = zeros(vidObj.Height,vidObj.Width,numFrames,'uint16');
iter = 1;
while hasFrame(vidObj)
    outimage(:,:,iter) = readFrame(vidObj);
    iter = iter+1;
end

mytiffwrite(outtifname,outimage)
% delete raw file, -f forces to delete write protected files
%     unix(sprintf('rm -f %s',tiffile))
if delete_file
    unix(sprintf('rm -f %s',videofile))
end

%%
end

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
