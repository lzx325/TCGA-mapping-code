if [ "$#" -lt 4 ]; then
    echo "Insufficient arguments" >&2
    exit 2
fi
workDir="${1%/}"
script_path="$workDir/anz/code/06.buildWhippletIndex/06.buildWhippletIndex.R"
cancer_type="$2"
whipplet_src="$3"
index="$4"
Rscript "$script_path" -w "$workDir" -c "$cancer_type" -d "$whipplet_src" -i "$index"