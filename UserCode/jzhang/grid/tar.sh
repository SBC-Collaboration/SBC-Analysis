#tar -C /bluearc/storage/SBC-17-data/$1 -czf /bluearc/storage/SBC-17-tgz/$1.tgz *
#tar -C /bluearc/storage/SBC-17-data/$1 -czf /bluearc/storage/SBC-17-tgz/$1.tgz *
#pushd /nashome/j/jzhang/Matlab/COUPPcode/Tools/$1
pushd /bluearc/storage/SBC-17-data/$1
tar -czf /bluearc/storage/SBC-17-tgz/$1.tgz *
