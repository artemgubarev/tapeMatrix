FILES=("err.log" "out.log" "tapeMatrix" "gen_err.log" "gen_out.log" "solution.txt" "msolution.txt")

for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        rm "$file"
        echo "�����: $file"
    else
        echo "���� �� ������: $file"
    fi
done