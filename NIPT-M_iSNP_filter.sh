#!/bin/bash

# Set input GCVF file
GVCF='/dssg/home/acct-medkwf/medkwf4/results/NIPT-M/merge/combined.raw.all.vcf'

# Set output file
output='/dssg/home/acct-medkwf/medkwf4/results/NIPT-M/merge/filter.vcf'
rm -fr ${output}
touch ${output}

# Function making float division calculation
calc(){ awk "BEGIN { print "$*" }"; }

IFS=$'\n'
count_row=0
FetalSum=0
for line in `< ${GVCF}`
do
  parental=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);print a[1]}')
  maternal=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);print a[1]}')
  fetal=$(echo $line | tr -d '\n' | awk '{n=split($12,a,/:/);print a[1]}')
  if [ "$parental" = "0/0" ] && [ "$maternal" = "1/1" ] && [ "$fetal" = "0/1" ]; then
    echo $line
    WTreads=$(echo $line | tr -d '\n' | awk '{n=split($12,a,/:/);split(a[2],b,/,/);print b[1]}')
    SNPreads=$(echo $line | tr -d '\n' | awk '{n=split($12,a,/:/);split(a[2],b,/,/);print b[2]}')
    depth=`calc $SNPreads+$WTreads`
    AF=`calc 1-$SNPreads/$depth`
    FetalFraction=`calc 2*$AF`
    echo -e "Fetal fraction is $FetalFraction \n"
    FetalSum=`calc $FetalFraction+$FetalSum`
    count_row=`calc $count_row+1`
    echo $line >> ${output}
  elif [ "$parental" = "1/1" ] && [ "$maternal" = "0/0" ] && [ "$fetal" = "0/1" ]; then
    echo $line
    WTreads=$(echo $line | tr -d '\n' | awk '{n=split($12,a,/:/);split(a[2],b,/,/);print b[1]}')
    SNPreads=$(echo $line | tr -d '\n' | awk '{n=split($12,a,/:/);split(a[2],b,/,/);print b[2]}')
    depth=`calc $SNPreads+$WTreads`
    AF=`calc $SNPreads/$depth`
    FetalFraction=`calc 2*$AF`
    echo -e "Fetal fraction is $FetalFraction \n"
    FetalSum=`calc $FetalFraction+$FetalSum`
    count_row=`calc $count_row+1`
    echo $line >> ${output}
  fi
done

FetalFraction_average=`calc $FetalSum/$count_row`
echo -e "Summary:\nTotal rows: ${count_row}\nFetal fraction: ${FetalFraction_average}"
