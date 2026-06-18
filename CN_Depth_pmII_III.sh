##########################################################################################
# PM II/III copy number assessment from depth of coverage

# Generate per-base coverage across PM II/III region
for file in *.bam
do
  bedtools genomecov -ibam "$file" -d |
    awk '$1 == "Pf3D7_14_v3" && $2 > 260000 && $2 < 327000' \
    > "${file%.bam}.PMs.coverage"
done

# Create output with sample names and copy number estimate
echo -e "sample\tPM23_depth\tflank_depth\tPM23_copy_number" > PMs.copy_number.tsv

for file in *.PMs.coverage
do
  sample=${file%.PMs.coverage}

  pm_depth=$(awk '$2>=289570 && $2<=298796 {tot+=$3; cnt++}
                  END {if (cnt>0) print tot/cnt; else print "NA"}' "$file")

  flank_depth=$(awk '$2<289570 || $2>298796 {tot+=$3; cnt++}
                     END {if (cnt>0) print tot/cnt; else print "NA"}' "$file")

  copy_number=$(awk -v pm="$pm_depth" -v flank="$flank_depth" \
                    'BEGIN {if (flank>0) print pm/flank; else print "NA"}')

  echo "$pm_depth" >> 1.coverage
  echo "$flank_depth" >> 2.coverage
  echo -e "${sample}\t${pm_depth}\t${flank_depth}\t${copy_number}" >> PMs.copy_number.tsv
done

# Keep your original two-column output too
paste 1.coverage 2.coverage > PMs.coverage


cat PMs.copy_number.tsv
