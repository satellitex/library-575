#!/bin/bash
dir="FinalZukky"
nkf -w --overwrite ./$dir/*

cd ./$dir
for f in *; do
    echo "//================= $f =====================" > ../result/$f
    cat $f >> ../result/$f
    echo "
    " >> ../result/$f
done

cd ../
fn="cat "
touch $fn
for f in ./result/*; do
    echo $f
    fn="$fn $f"
done

$fn > print_result.txt

