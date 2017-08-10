% read original images 
myfold = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/comp1/'
myfiles = dir([myfold,'*.h5'])
outfolder = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/compressed_data_for_training/'
mkdir(outfolder)
for ii=3:2:length(myfiles)
    myfile = fullfile(myfold,myfiles(ii).name);
    hinfo = hdf5info(myfile);
    imgdata = permute(hdf5read(myfile,hinfo.GroupHierarchy.Datasets.Name),[3 2 4 1]);
    dims = size(imgdata);
    
    iter = 0;
    Comp = [5 10 20 40 60 80]
    etime = zeros(1,18);
    parfor iter =1:6
        comp = Comp(iter);
%         tic
        videoname = fullfile(outfolder,sprintf('%s%d.%s',myfiles(ii).name(1:end-4),comp,'mj2'));
        v = VideoWriter(videoname,'Motion JPEG 2000');
        v.CompressionRatio = comp;
        v.MJ2BitDepth = 16;
        open(v)
%         writeVideo(v,reshape(imgdata,[size(imgdata,1),size(imgdata,2),1,size(imgdata,3)]))
        writeVideo(v,reshape(imgdata,[size(imgdata,1),size(imgdata,2),1,size(imgdata,3)]))
        close(v)
%         etime((iter-1)*3+1) = toc;
        
        %% read it back
%         tic
        v = VideoReader(videoname);
        video = zeros(size(imgdata),'uint16');
        k=1;
        while hasFrame(v)
            video(:,:,k) = readFrame(v);
            k = k+1;
        end
%         etime((iter-1)*3+2) = toc;
        %%
        clc
%         tic
        % write it as chuncked hdf5
        outfile = fullfile(outfolder,sprintf('%s%d.%s',myfiles(ii).name(1:end-4),comp,'h5'))
        
        video_out = permute(video,[2 1 3]);
        dims = size(video_out);
        video_out = reshape(video_out,[1 dims]);
        dims = [size(video_out)];
        h5create(outfile,'/exported_data',dims,'Datatype','uint16','ChunkSize',[1,64,63,63],'Deflate',1)
        h5write(outfile,'/exported_data', video_out)
%         etime((iter-1)*3+3) = toc;
    end
%     Etime(ii) = etime;
end
%%
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
%find number of random characters to choose from
numRands = length(s); 
%specify length of random string to generate
sLength = 10;
infolder = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/compressed_data_for_training/'
logfolder = '/nobackup/mousebrainmicro/cluster/logs/'
filenames = {'3727-0','5999-0','7733-0'}
ilastikloc = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/ilastik-1.1.8.post1-Linux/run_ilastik.sh'

myfile = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/compressed_data_for_training/compressiontest.sh'
% test different classifiers
fid = fopen(myfile,'w');
numcores = 8;
memsize = numcores*7.5*1000;

for classifier = [5 10 20 40 60 80]
    ilpfile = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/compressed_data_for_training/compressionTest_cmp%02d.ilp',classifier);
    
    for ii=1:3
        %generate random string
        randString = s( ceil(rand(1,sLength)*numRands) );
        name = sprintf('ilp_%05d-cmp_%d-%s',ii,classifier,randString);
        logfile=fullfile(logfolder,sprintf('ilastikclassifierrun_%05d-%s.txt',ii,randString));
        infiles=fullfile(infolder,sprintf('%s_comp_%02d.h5',filenames{ii},classifier));
        outputformat=sprintf('''%s%s_comp_%02d_Probabilities.h5''',infolder,filenames{ii},classifier);
        args = sprintf('''%s --headless --logfile=%s --project=%s --output_format=hdf5 --output_filename_format=%s %s''',...
            ilastikloc,logfile,ilpfile,outputformat,infiles);
        mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d qsub -pe batch %d -l short=true -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,memsize,numcores,name,args);
        fwrite(fid,mysub);
    end
end
unix(sprintf('chmod +x %s',myfile));
%%
% read classifier result
outloc = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/compressed_data_for_training/'
IP = []
iter =0;
for classifier = [1 5 10 20 40 60 80]
    iter = iter+1;
    for ii=1:3
        outputformat=sprintf('%s%s_comp_%02d_Probabilities.h5',outloc,filenames{ii},classifier);
        h5info = hdf5info(outputformat);
        IP{iter}{ii} = hdf5read(outputformat,h5info.GroupHierarchy.Datasets.Name);
    end
end
%%
% close all
iter = 0;
mm=0;
for classifier = [1 5 10 20 40 60 80]
    iter = iter+1;
    for ii=1:3
        mm = mm+1;
        It = squeeze(IP{iter}{ii}(1,:,:,:));
        Im1 = squeeze(max(It,[],3));
        Im2 = squeeze(max(It.*(It>.5),[],3));
        figure(100+9)
        subplot(3,7,iter+(ii-1)*7)
        imshow(Im1',[])
        title(sprintf('CL: %d, IM: %d',classifier,ii))
        figure(200+9)
        subplot(3,7,iter+(ii-1)*7)
        imshow(Im2',[])
        title(sprintf('CL: %d, IM: %d',classifier,ii))
        drawnow
    end
end
%%
iter = 0;
mm=0;
for classifier = [1 5 10 20 40 60 80]
    iter = iter+1;
    for ii=1:3
        mm = mm+1;
        It = squeeze(IP{iter}{ii}(1,:,:,:));
        It(:,:,1:10) = 0;
        
        Im1 = squeeze(max(It,[],3));
        
        CC = bwconncomp(It>.5);
        lengths = cellfun(@length,CC.PixelIdxList);
        M = zeros(size(It));
        for jj = 1:CC.NumObjects
            if lengths(jj)>1e3
                M(CC.PixelIdxList{jj})=1;
            end
        end
        
        Im2 = squeeze(max(It.*M,[],3));
        figure(100+9)
        subplot(3,7,iter+(ii-1)*7)
        imshow(Im1',[])
        title(sprintf('CL: %d, IM: %d',classifier,ii))
        figure(200+9)
        subplot(3,7,iter+(ii-1)*7)
        imshow(Im2',[])
        title(sprintf('CL: %d, IM: %d',classifier,ii))
        drawnow
    end
end
%%

iter = 0;
mm=0;
for classifier = [1 5 10 20 40 60 80]
    iter = iter+1;
    for ii=1:3
        mm = mm+1;
        It = squeeze(IP{iter}{ii}(1,:,:,:));
        It(:,:,1:10) = 0;
        CC = bwconncomp(It>.5);
        lengths = cellfun(@length,CC.PixelIdxList);
        M = zeros(size(It));
        for jj = 1:CC.NumObjects
            if lengths(jj)>1e3
                M(CC.PixelIdxList{jj})=1;
            end
        end
        
        Im2 = squeeze(max(It.*M,[],3));
        IO{iter,ii} = Im2;
    end
end
%%
tic
v = VideoReader('/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/compressed_data_for_training/3727-0_comp_80.mj2');
video = zeros(size(imgdata),'uint16');
k=1;
while hasFrame(v)
    video(:,:,k) = readFrame(v);
    k = k+1;
end
toc
%%
























