if [ "$#" -lt 4 ]; then
    echo "Insufficient arguments" >&2
    exit 2
fi
workDir="${1%/}"
script_path="$workDir/anz/code/07.WhippletQuant/07.WhippletQuant.R"
cancer_type="$2"
WhippletQuantjls="$3"
index="$4"
Rscript "$script_path" -w "$workDir" -c "$cancer_type" --WhippletQuantjls "$WhippletQuantjls"  --index "$index"