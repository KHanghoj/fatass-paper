awk '{print 1/NF, 1/NF, 1/NF}' ../res_inclafr/admix/K3_s1/sub0.q > flat.q
ln -f -s ../res_inclafr/admix/K3_s1/sub*.npy .
ln -f -s ../res_inclafr/haplonet_split/allchrom.sub*.npy .
HAPLONET=/home/jonas/miniconda3/bin/haplonet

mkdir -p res

threads=40

for x in 0 2 4
do
	q=flat.q
	f=sub${x}.f.npy
	l=allchrom.sub${x}.loglike.npy
	out=res/${x}_flatq_1000
	${HAPLONET} fatash --no_optim --alpha_save --alpha 1000 --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
done

for x in 0 2 4
do
	q=sub${x}.q
	f=sub${x}.f.npy
	l=allchrom.sub${x}.loglike.npy
	out=res/${x}_realq_1000
	${HAPLONET} fatash --no_optim --alpha_save --alpha 1000 --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
done

for x in 0 2 4
do
	q=sub${x}.q
	f=sub${x}.f.npy
	l=allchrom.sub${x}.loglike.npy
	out=res/${x}_realq_0.01
	${HAPLONET} fatash --no_optim --alpha_save --alpha 0.01 --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
done

for x in 0 2 4
do
	q=sub${x}.q
	f=sub${x}.f.npy
	l=allchrom.sub${x}.loglike.npy
	out=res/${x}_realq_scaled
	if [[ x -eq 0 ]]
	then 
		alpha=0.01
	else
		alpha=$( echo  "scale=6; 0.01 / $x" | bc)
	fi
	echo $alpha
	${HAPLONET} fatash --no_optim --alpha_save --alpha $alpha --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
done

for x in 0 2 4
do
	q=sub${x}.q
	f=sub${x}.f.npy
	l=allchrom.sub${x}.loglike.npy
	out=res/${x}_realq_est
	${HAPLONET} fatash --alpha_save --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
done
