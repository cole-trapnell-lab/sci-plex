# awk variable PCR_COMBO needs to be defined
# input 1: sci-RNA-seq.RT.oligos
# input 2: combinatorial.indexing.key
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
	rt_index[$2] = $3;
	rt_name[$1] = $2;
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
        # special case: a file that's just a sample name and nothing else
        if ($0 ~ /^#/) next;
        if (NF == 1) {
            single_sample = $1;
        } else if (p5_row >= $2 && p5_row <= $3 && p5_col >= $4 && p5_col <= $5) {
            rt_row_min[$1] = $6;
            rt_row_max[$1] = $7;
            rt_col_min[$1] = $8;
            rt_col_max[$1] = $9;
        }
    } else {
        read_num++;
	header=$0
        getline;
        rt_barcode = substr($1, 1, 10);
        polyA = substr($1, 12, 4);
	#printf "%s\n", $0;
        if (rt_barcode in rt_well && polyA == "AAAA") {
            hits++;

            this_rt_well = rt_well[rt_barcode];
	    this_rt_index = rt_index[rt_name[this_rt_well]];
            rt_row = substr(this_rt_well, 1, 1);
            rt_col = substr(this_rt_well, 2, 2);

            this_sample = "?";

            if (single_sample != "") {
                this_sample = single_sample;
            } else {
                for (sample in rt_col_min) {
                    if ((rt_row_min[sample] == "NA" || rt_row >= rt_row_min[sample]) &&
                        (rt_row_max[sample] == "NA" || rt_row <= rt_row_max[sample]) &&
                        (rt_col_min[sample] == "NA" || rt_col >= rt_col_min[sample]) &&
                        (rt_col_max[sample] == "NA" || rt_col <= rt_col_max[sample])) {
                        this_sample = sample;
                    }
                }
            }

            printf "%s|%s|%s\n", header, this_rt_well, this_rt_index;
            getline;
            getline;
        } else {
            getline;
            getline;
        }
    }
} END {
}

