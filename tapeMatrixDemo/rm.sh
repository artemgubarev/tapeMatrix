FILES=("err.log" "out.log" "tapeMatrix" "gen_err.log" "gen_out.log" "solution.txt" "msolution.txt")

for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        rm "$file"
        echo "”далЄн: $file"
    else
        echo "‘айл не найден: $file"
    fi
done