# awk variable PCR_COMBO needs to be defined
# input 1: sci-RNA-seq.RT.oligos
# input 2: sci-RNA-seq.RT.samples
# input 3: paste applied to read 1 and read 2 files (gunzipped in subshells)

BEGIN {
    read_num = 0;
    hits = 0;

    bases[1] = "A";
    bases[2] = "C";
    bases[3] = "G";
    bases[4] = "T";

    p5_row = substr(PCR_COMBO, 1, 1);
    p5_col = substr(PCR_COMBO, 2, 2);
    p7_row = substr(PCR_COMBO, 4, 1);
    p7_col = substr(PCR_COMBO, 5, 2);

    single_sample = "";
} {
    if (ARGIND == 1) {
        rt_well[$2] = $1;
        for (i = 1; i <= length($2); i++) {
            for (j = 1; j <= 4; j++) {
                mismatch = "";
                if (i > 1)
                    mismatch = mismatch substr($2, 1, i-1);
                mismatch = mismatch bases[j];
                if (i < length($2))
                    mismatch = mismatch substr($2, i+1, length($2)-i);
                if (!(mismatch in rt_well))
                    rt_well[mismatch] = $1;
            }
        }
    } else if (ARGIND == 2) {
        if ($0 ~ /^#/) next;
        if (FNR == 1 && NF == 1) {
            # special case: a file that's just a sample name and nothing else
            single_sample = $1;
            nextfile;
        } else {
            RT_sample[$1] = $2;
        }
    } else {
        read_num++;
        getline;
        rt_barcode = substr($1, 9, 10);
        umi = substr($1, 1, 8);

        if (rt_barcode in rt_well) {
            hits++;

            this_rt_well = rt_well[rt_barcode];
            this_sample = "?";

            if (single_sample != "") {
                this_sample = single_sample;
            } else if (this_rt_well in RT_sample) {
                this_sample = RT_sample[this_rt_well];
            }

            printf "@%s_%s|%s|%s|%s|%s|%s\n",
                PCR_COMBO, read_num, this_sample,
                substr(PCR_COMBO, 1, 3), substr(PCR_COMBO, 4, 3),
                this_rt_well, umi;
            printf "%s\n", $2;
            getline;
            printf "+\n";
            getline;
            printf "%s\n", $2;
        } else {
            getline;
            getline;
        }
    }
} END {
    printf("%d\t%d\t%.3f\t(RT barcode matches, total reads, proportion)\n",
        hits, read_num, hits / read_num) > "/dev/stderr";
}

