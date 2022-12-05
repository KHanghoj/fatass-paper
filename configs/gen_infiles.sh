grep -w -f na_inclafr.popnames integrated_call_samples_v2.20130502.ALL.panel | awk '{print $1}' > na_inclafr.samples
grep -w -f na_inclafr.popnames integrated_call_samples_v2.20130502.ALL.panel | awk '{print $2}' > na_inclafr.pops
grep -w -f na_inclafr_eas.popnames integrated_call_samples_v2.20130502.ALL.panel | awk '{print $1}' > na_inclafr_eas.samples
grep -w -f na_inclafr_eas.popnames integrated_call_samples_v2.20130502.ALL.panel | awk '{print $2}' > na_inclafr_eas.pops
