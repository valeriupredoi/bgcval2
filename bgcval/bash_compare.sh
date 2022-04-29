module load jaspy/2.7
x=0
#while [ 1 -eq 1 ];
while [ $x -le 35 ]
do
    ./analysis_compare.py
    rsync -av /home/users/ldemora/workspace/ukesm-validation/CompareReports2 /gws/nopw/j04/esmeval/public/CompareReports/v2/.
    echo "sleeping" $x
    x=$(( $x + 1 ))
    #leep 60000;
    sleep $((($(date -f - +%s- <<<$'06:00 tomorrow\nnow')0)%86400))
done
