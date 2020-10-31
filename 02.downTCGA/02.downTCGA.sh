if [ "$#" -lt 2 ]; then
    echo "Insufficient arguments" >&2
    exit 2
fi
workDir="${1%/}"
cancer="$2"
script_path="$workDir/anz/code/02.downTCGA/02.downTCGA.R"
Rscript "$script_path" -w "$workDir" -c "$cancer" -t 10 
