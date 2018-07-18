if [[ -f /grid/fermiapp/products/common/etc/setups.sh ]]; then
	source /grid/fermiapp/products/common/etc/setups.sh
	setup ifdhc
fi

code_dir=/bluearc/storage/recon/devel/code
pnfs_dir=/pnfs/coupp/persistent/runlists/SBC-17

if [ -f $code_dir/SBCcode.tar ]; then
	rm $code_dir/SBCcode.tar
fi

cd $code_dir/
tar -cvf $code_dir/SBCcode.tar SBCcode \
    --exclude "UserCode/jzhang/log" \
    --exclude-vcs

if ifdh ls $pnfs_dir/SBCcode.tar; then
    echo "ifdh rm $pnfs_dir/SBCcode.tar"
	ifdh rm $pnfs_dir/SBCcode.tar
fi
echo "copying SBCcode.tar to /pnfs"
ifdh cp $code_dir/SBCcode.tar \
    $pnfs_dir/SBCcode.tar

