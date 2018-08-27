if [ -f /coupp/data/home/coupp/sbc_grid_test/SBCcode.tar ]; then
	rm /coupp/data/home/coupp/sbc_grid_test/SBCcode.tar
fi
cd /bluearc/storage/recon/devel/code/
tar -cvf /coupp/data/home/coupp/sbc_grid_test/SBCcode.tar SBCcode \
    --exclude "UserCode/jzhang/log" \
    --exclude-vcs

if ifdh ls /pnfs/coupp/persistent/runlists/SBC-17/SBCcode.tar; then
    echo "ifdh rm /pnfs/coupp/persistent/runlists/SBC-17/SBCcode.tar"
	ifdh rm /pnfs/coupp/persistent/runlists/SBC-17/SBCcode.tar
fi
echo "copying SBCcode.tar to /pnfs"
ifdh cp /coupp/data/home/coupp/sbc_grid_test/SBCcode.tar \
    /pnfs/coupp/persistent/runlists/SBC-17/SBCcode.tar

