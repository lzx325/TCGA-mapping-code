if [ "$#" -lt 3 ]; then
    echo "Insufficient arguments" >&2
    exit 2
fi
workDir="${1%/}"
script_path="$workDir/anz/code/05.STAR_2ndAlign/05.STAR_2ndAlign.R"
Rscript "$script_path" -w "$workDir" -c "$2" -i "$3" --thread "${4:-30}"
