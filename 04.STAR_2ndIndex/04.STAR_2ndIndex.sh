if [ "$#" -lt 2 ]; then
    echo "Insufficient arguments" >&2
    exit 2
fi
workDir="${1%/}"
script_path="$workDir/anz/code/04.STAR_2ndIndex/04.STAR_2ndIndex.R"
Rscript "$script_path" -w "$workDir" -c "$2" --thread "${3:-30}"