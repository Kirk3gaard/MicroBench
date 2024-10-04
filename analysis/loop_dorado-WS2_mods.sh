#!/bin/bash

# You need to stop minknow and other stuff using the GPUs
#sudo service minknow stop
sudo service doradod stop

# Input
FLOWCELLID=PAW77640;
INPUTDIR="/data/2024-10-01_np_PAW77640/no_sample_id/20241001_1131_P2S-00692-A_PAW77640_9c9e38df/pod5/"

############################
DORADOVERSION="0.7.3" # data model version 5.0.0
BASECALLINGMODEL="sup" # Choose between fast, hac, sup
MODELVERSION="5.0.0"
DORADO="/data/software/dorado-"$DORADOVERSION"-linux-x64/bin/dorado"

MODS="4mC_5mC 6mA" # https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models
#

#dorado basecaller --device "cuda:0" $BASECALLINGMODEL $INPUTDIR --modified-bases $MODS > $OUTPUTFILE.mod.bam

mkdir -p current_file/
   
for f in /data/2024-10-01_np_PAW77640/no_sample_id/20241001_1131_P2S-00692-A_PAW77640_9c9e38df/pod5/*.pod5; do
   FILEID=$(echo $f | sed -E 's/.*_cddc0ba5_([0-9]+).pod5/\1/')
   ln -s $f current_file/this.pod5
   OUTPUTFILE=$FLOWCELLID".dorado"$DORADOVERSION".bm"$MODELVERSION"."$BASECALLINGMODEL.mod4mC_5mC_6mA.$FILEID.bam
   dorado basecaller --device "cuda:0" $BASECALLINGMODEL current_file/ --modified-bases $MODS > $OUTPUTFILE
   rm current_file/this.pod5
done



