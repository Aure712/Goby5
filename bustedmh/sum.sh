cat Relax_T_only/* | awk -F'\t' '$9 < 0.05' | grep 'Success' > sumT_p0.05.txt
