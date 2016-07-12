function deployment
%qsub -pe batch 4 -l short=true -N tile_test -j y -o ~/logs -b y -cwd -V './compiledfiles_mytest/mytest > output_mytest.log'
%%
% LAZYFLOW_THREADS=8 LAZYFLOW_TOTAL_RAM_MB=60000 ilastik_cluster --headless --logfile=/nobackup/mousebrainmicro/cluster/2015-06-19-GN1-Takako/logs/ilastikclassifierrun.txt --project=/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/GN1.ilp --output_format=hdf5 --output_filename_format='/nobackup/mousebrainmicro/cluster/2015-06-19-GN1-Takako/classifier_output/{nickname}_Probabilities.h5' /nobackup/mousebrainmicro/cluster/2015-06-19-GN1-Takako/ch1_h5/IM-0001_idx-03403_ch-1.h5


functionname = ''
numcores = 8;
memsize = numcores*7.5*1000;

logfolder = '/nobackup/mousebrainmicro/cluster/2015-06-19-GN3-Takako/logs/'
mkdir(logfolder)

out = '/tier2/mousebrainmicro/mousebrainmicro/cluster/compressionExperiment/out/'
mkdir(out)

myfile = 'cluster_ilastik_GN3_compexps.sh';

s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
%find number of random characters to choose from
numRands = length(s);
%specify length of random string to generate
sLength = 10;

ilastikloc = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/ilastik-1.1.8.post1-Linux/run_ilastik.sh'
%[2648    3392    3752    5666    7362    8368    8369    8716    8717]
fid = fopen(myfile,'w');
idxTiles = [2648    3392    3752    5666    7362    8368    8369    8716    8717];
% for ii=1:length(myfiles)
for comp = [1 5 10 20 40 80 160];
    outfolder_ = sprintf('/tier2/mousebrainmicro/mousebrainmicro/cluster/compressionExperiment/comp%02d/mj2',comp);
    % outfolder_ = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/nobackup_tmp/comp_%d/GN3',comp);
    % read an image
    for idx = 1:length(idxTiles);
        %%
        [comp idx]
        outfolder = outfolder_;
        for ch = [0 1]
            outfolder = sprintf('/tier2/mousebrainmicro/mousebrainmicro/cluster/compressionExperiment/comp%02d/h5/',comp);
            myfiles = fullfile(outfolder,sprintf('IM-%04d_idx-%05d_ch-%d_comp_%d.h5',idx,idxTiles(idx),ch,comp));
            
            %generate random string
            if comp==1
                ilpfile = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/GN1.ilp';
            else
                ilpfile = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/software/training/GN1-training/data/compressed_data_for_training/compressionTest_cmp%02d.ilp',comp);
            end
            randString = s( ceil(rand(1,sLength)*numRands) );
            name = sprintf('ilp_%05d-%s',idx,randString);
            logfile=fullfile(logfolder,sprintf('ilastikclassifierrun_%05d_%d-%s.txt',idx,idx,randString));
            infiles=myfiles;%
            outputformat=sprintf('''%s_Probabilities.h5''',fullfile(out,myfiles(1:end-3)));
            args = sprintf('''%s --headless --logfile=%s --project=%s --output_format=hdf5 --output_filename_format=%s %s''',...
                ilastikloc,logfile,ilpfile,outputformat,infiles);
            
            mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d qsub -pe batch %d -l short=true -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,memsize,numcores,name,args);
            fwrite(fid,mysub);
        end
        %%
        %     mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d qsub -pe batch %d -l short=true -N %s -j y -o ~/logs -b y -cwd -V %s\n',numcores,memsize,numcores,name,args);
        %     outputformat=sprintf('''/nobackup/mousebrainmicro/cluster/2015-06-19-GN1-Takako/classifier_output/{nickname}_Probabilities.h5''');
        %     args = sprintf('''LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d %s --headless --logfile=%s --project=%s --output_format=hdf5 --output_filename_format=%s %s''',...
        %         numcores,memsize,ilastikloc,logfile,ilpfile,outputformat,infiles);
        
    end
end
        unix(sprintf('chmod +x %s',myfile));

end
