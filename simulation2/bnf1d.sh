#!/bin/bash

dirname='dir-dat'
#command='./d2'
command='./d1'
#p='d2'				# prefix
p='d1'				# prefix
a=0.01				# alpha: learning rate
k=3				# k-near bit
m=1				# matrix scale
s=0				# random seed
r=1				# rayout   1:multi_image, 2:2d_data, 6:2d_all
t=1				# terminal 0:eps, 1:gif, 2:x11

ids=`cat <<EOF
0
10
20
30
40
50
60
70
80
90
100
1000
2000
3000
4000
5000
6000
7000
8000
9000
10000
9900
9910
9920
9930
9940
9950
9960
9970
9980
9990
EOF`

g2e(){
    files=`ls ${1}/*.gif`
    for file in files; do
	for id in ids; do
	    if [[ "${file}" =~ t${id}\.gif$ ]]; then
		echo "${file} is converting gif to eps."
		convert -crop 639x129+0+170 ${file} ${file%gif}eps
	    fi
	done
    done
}

#for a in 0.01 0.1 0.5; do
#for k in 3 4 5; do
#    for m in 0.01 0.1 0.5; do
for s in 0 1 2 3 4; do
#	for r in 1 2 ; do 
    cmdopt=`printf '-a%4.2f -k%d -m%4.2f -s%d -r%d' ${a} ${k} ${m} ${s} ${r}`
    echo "cmdopt=${cmdopt}"
    #cmdopt="-a${a} -k${k} -m${m} -s${s} -r${r}"

    prefix=`echo ${cmdopt} | sed -e 's/[ -]//g' | sed -e 's/\./_/g'`
    prefix=`echo "${p}-${prefix}" | sed -e 's/\./_/g'`
    echo "prefix=${prefix}"
    odir="${dirname}/${prefix}"
    echo mkdir -p ${odir}
    echo ${cmdopt}
    echo ${command} ${cmdopt} -t${t}
#    echo convert `seq -f "${dirname}/${prefix}-%g.gif" 10 10 10000` "${dirname}/${prefix}-anime.gif"
    echo convert 'seq -f "${dirname}/${prefix}-%g.gif" 10 10 10000' "${dirname}/${prefix}-anime.gif"
    echo mv ${dirname}/${prefix}*.gif ${odir}
    echo g2e ${odir}
done

#done
#done
#done
