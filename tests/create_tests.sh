from="/home/filiper/Projetos/CpH_project/CpHMD-old/pentapeptides/NTR/control/OLD/pH"
pH="07.2"
to="penta_control/NTR"
base=`pwd`

rm -f ${to}/old_mocc
for rep in 1 2 3
do
    target=${from}${pH}_r${rep}
    
    cp $target/*.gro $to

    cd $to
    for i in `ls protein*.gro`
    do 
        grep -v SOL $i > ${i}_
        sed -i "s/2951/35/" ${i}_  #"s/1827/48/" ${i}_ 
        mv ${i}_ ${rep}_$i
        rm -f ${i}
    done
    rm -f *protein*50.gro

    for i in `ls $target/*.mocc`
    do 
        head -n 1 $i >> old_mocc
    done
    cd $base
done

#cat tmp.out | histog -r -0.05,0.05 -n 10
#sed -n '1~4p' tmp.out | histog -r -0.05,0.05 -n 10
#sed -n '2~4p' tmp.out | histog -r -0.05,0.05 -n 10
#sed -n '3~4p' tmp.out | histog -r -0.05,0.05 -n 10
#sed -n '4~4p' tmp.out | histog -r -0.05,0.05 -n 10