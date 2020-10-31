#!/usr/bin/env bash
set -e
cancers=("TCGA-CESC" "TCGA-UCEC")
server="anzheng@proxy.y1zhou.com"
local_root="/ibex/scratch/projects/c2066/anz/result"
remote_root="/mnt/csbl5/anzheng/alternative_splicing/kaust"
. ~/anz_config.sh
for c in "${cancers[@]}"
do
    echo "$c"
    ssh -p 27255 "$server"  "mkdir -p $remote_root/$c $remote_root/$c/03.STAR_1stAlign $remote_root/$c/06.buildWhippletIndex $remote_root/$c/07.WhippletQuant"
    sync_us "$local_root/$c/03.STAR_1stAlign/" "$remote_root/$c/03.STAR_1stAlign" 
    sync_us "$local_root/$c/06.buildWhippletIndex/" "$remote_root/$c/06.buildWhippletIndex" "--include=*.exons.tab.gz --exclude=* "
    sync_us "$local_root/$c/07.WhippletQuant/" "$remote_root/$c/07.WhippletQuant" 
done
