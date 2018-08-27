data_dir="/bluearc/storage/SBC-17-data"
cd $data_dir

RUN_LIST="/pnfs/coupp/persistent/runlists/SBC-17/list_1499442343"
#runs=`cat $RUN_LIST`
#for run in $runs; do
#    echo $run
#done

# LEN=`cat $RUN_LIST | wc -l`
# for LINE in `seq 69 $LEN`; do
#     run=`head -$LINE $RUN_LIST | tail -1`
#     pushd $run
#     md5sum *.* > /nashome/j/jzhang/SBCcode/UserCode/jzhang/log/$run.md5
#     dirlist=`ls -d */ | sort -V`
#     for d in $dirlist; do
#         md5sum $d* >> /nashome/j/jzhang/SBCcode/UserCode/jzhang/log/$run.md5
#     done
#     popd
# done

##############################
for run in 20170622_0 20170623_3 20170624_4 20170625_0 20170627_0 20170629_0 20170630_1 20170630_8 20170702_4 20170703_7 20170704_3 20170704_4
do
    pushd $run
    md5sum *.* > /nashome/j/jzhang/SBCcode/UserCode/jzhang/log/$run.md5
    dirlist=`ls -d */ | sort -V`
    for d in $dirlist; do
        md5sum $d* >> /nashome/j/jzhang/SBCcode/UserCode/jzhang/log/$run.md5
    done
    popd
done


