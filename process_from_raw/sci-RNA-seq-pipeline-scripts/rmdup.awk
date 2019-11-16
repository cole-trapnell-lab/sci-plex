BEGIN {
    prev_pos = 0;
    n_prev_at_pos = 0;
    prev_rt_well[1] = "";
    prev_umi[1] = "";
} {
    if (substr($0, 1, 1) == "@") {
        print; next;
    }

    split($1, arr, "|");

    pos = $4;
    rt_well = arr[5];
    umi = arr[6];

    if (pos == prev_pos) {
        for (i = 1; i <= n_prev_at_pos; i++) {
            #print i "\t" prev_rt_well[i] "\t" prev_umi[i];

            if (rt_well == prev_rt_well[i]) {
                n_mismatch = 0;
                for (j = 1; j <= 8; j++) {
                    if (substr(umi, j, 1) != substr(prev_umi[i], j, 1))
                        n_mismatch++;
                }

                if (n_mismatch == 0)
                    next;
                else if (n_mismatch == 1) {
                    n_prev_at_pos++;
                    prev_rt_well[n_prev_at_pos] = rt_well;
                    prev_umi[n_prev_at_pos] = umi;
                    next;
                }
            }
        }

        print;
        n_prev_at_pos++;
        prev_rt_well[n_prev_at_pos] = rt_well;
        prev_umi[n_prev_at_pos] = umi;
    } else {
        print;

        prev_pos = pos;
        n_prev_at_pos = 1;
        prev_rt_well[n_prev_at_pos] = rt_well;
        prev_umi[n_prev_at_pos] = umi;
    }
}

