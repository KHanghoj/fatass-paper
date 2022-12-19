f=integrated_call_samples_v2.20130502.ALL.panel
grep -w -f na_inclafr.popnames ${f} | awk '{print $1}' > na_inclafr.samples
grep -w -f na_inclafr.popnames ${f} | awk '{print $2}' > na_inclafr.pops
grep -w -f na_inclafr_eas.popnames ${f} | awk '{print $1}' > na_inclafr_eas.samples
grep -w -f na_inclafr_eas.popnames ${f} | awk '{print $2}' > na_inclafr_eas.pops
grep -w -f aa.popnames ${f} | grep -v -w -f excl_related.txt | awk '{print $1}' > aa.samples
grep -w -f aa.popnames ${f} | grep -v -w -f excl_related.txt | awk '{print $2}' > aa.pops
