FILES=("err.log" "out.log" 
        "tapeMatrix" 
        "gen_err.log" "gen_out.log" 
        "solution.txt" "msolution.txt"
        "mout.log" "merr.log")

for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        rm "$file"
        echo "delete: $file"
    else
        echo "not found: $file"
    fi
done