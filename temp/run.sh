d=../res_inclafr
awk '{print 1/NF, 1/NF, 1/NF}' ${d}/admix/K3_s1/sub0.q > flat.q
ln -f -s ${d}/admix/K3_s1/sub*.npy .
ln -f -s ${d}/haplonet_split/allchrom.sub*.npy .
# HAPLONET=/home/jonas/miniconda3/bin/haplonet
HAPLONET=/home/krishang/mambaforge/envs/haplonet/bin/haplonet
mkdir -p res

threads=40

# for x in 0 2 4
# do
# 	q=flat.q
# 	f=${d}/admix/K3_s1/sub${x}.f.npy
# 	l=${d}/haplonet_split/allchrom.sub${x}.loglike.npy
# 	out=res/${x}_flatq_1000
# 	${HAPLONET} fatash --no_optim --alpha_save --alpha 1000 --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
# done

# for x in 0 2 4
# do
# 	q=${d}/admix/K3_s1/sub${x}.q
# 	f=${d}/admix/K3_s1/sub${x}.f.npy
# 	l=${d}/haplonet_split/allchrom.sub${x}.loglike.npy
# 	out=res/${x}_realq_1000
# 	${HAPLONET} fatash --no_optim --alpha_save --alpha 1000 --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
# done

# for x in 0 2 4
# do
# 	q=${d}/admix/K3_s1/sub${x}.q
# 	f=${d}/admix/K3_s1/sub${x}.f.npy
# 	l=${d}/haplonet_split/allchrom.sub${x}.loglike.npy
# 	out=res/${x}_realq_0.01
# 	${HAPLONET} fatash --no_optim --alpha_save --alpha 0.01 --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
# done

# for x in 0 2 4
# do
# 	q=${d}/admix/K3_s1/sub${x}.q
# 	f=${d}/admix/K3_s1/sub${x}.f.npy
# 	l=${d}/haplonet_split/allchrom.sub${x}.loglike.npy
# 	out=res/${x}_realq_scaled
# 	if [[ x -eq 0 ]]
# 	then 
# 		alpha=0.01
# 	else
# 		alpha=$( echo  "scale=6; 0.01 / $x" | bc)
# 	fi
# 	echo $alpha
# 	${HAPLONET} fatash --no_optim --alpha_save --alpha $alpha --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
# done

# for x in 0 2 4
# do
# 	q=${d}/admix/K3_s1/sub${x}.q
# 	f=${d}/admix/K3_s1/sub${x}.f.npy
# 	l=${d}/haplonet_split/allchrom.sub${x}.loglike.npy
# 	out=res/${x}_realq_scaled0005
# 	if [[ x -eq 0 ]]
# 	then 
# 		alpha=0.005
# 	else
# 		alpha=$( echo  "scale=6; 0.005 / $x" | bc)
# 	fi
# 	echo $alpha
# 	${HAPLONET} fatash --no_optim --alpha_save --alpha $alpha --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
# done

# for x in 0 2 4
# do
# 	q=${d}/admix/K3_s1/sub${x}.q
# 	f=${d}/admix/K3_s1/sub${x}.f.npy
# 	l=${d}/haplonet_split/allchrom.sub${x}.loglike.npy
# 	out=res/${x}_realq_scaled0001
# 	if [[ x -eq 0 ]]
# 	then 
# 		alpha=0.001
# 	else
# 		alpha=$( echo  "scale=6; 0.001 / $x" | bc)
# 	fi
# 	echo $alpha
# 	${HAPLONET} fatash --no_optim --alpha_save --alpha $alpha --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
# done

# for x in 0 2 4
# do
# 	q=${d}/admix/K3_s1/sub${x}.q
# 	f=${d}/admix/K3_s1/sub${x}.f.npy
# 	l=${d}/haplonet_split/allchrom.sub${x}.loglike.npy
# 	out=res/${x}_realq_scaled00001
# 	if [[ x -eq 0 ]]
# 	then 
# 		alpha=0.0001
# 	else
# 		alpha=$( echo  "scale=8; 0.0001 / $x" | bc)
# 	fi
# 	echo $alpha
# 	${HAPLONET} fatash --no_optim --alpha_save --alpha $alpha --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
# done

for x in 0 2 4
do
	q=${d}/admix/K3_s1/sub${x}.q
	f=${d}/admix/K3_s1/sub${x}.f.npy
	l=${d}/haplonet_split/allchrom.sub${x}.loglike.npy
	bp=${d}/basepos/allchrom_sub${x}.positions.txt
	out=res/${x}_realq_0.0001windows
	alpha=0.0001
	${HAPLONET} fatash --no_optim --window_save --alpha_save --alpha $alpha --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads} -w <( cut -f 2 ${bp})
done

# for x in 0 2 4
# do
# 	q=${d}/admix/K3_s1/sub${x}.q
# 	f=${d}/admix/K3_s1/sub${x}.f.npy
# 	l=${d}/haplonet_split/allchrom.sub${x}.loglike.npy
# 	out=res/${x}_realq_est
# 	${HAPLONET} fatash --alpha_save --like ${l} --prop ${q} --freq ${f}  --out ${out} --threads ${threads}
# done
