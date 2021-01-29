for d in rejoin merge scripts; do
    pushd $d && make clean && make && popd
done
