#!/bin/bash
path_name=$(pwd)
#find ./ -name "*.gz" -exec basename {} \; >1.txt #只列出最后的名字
find ./ -name "*.gz" >1.txt
for i in `cat 1.txt`;
do
#  echo ${i%%.*}.fg.gz
#  echo ${i%%.*}
  temp=${i#*/}  #删除从左边开始的第一个/后面的
#  echo $temp
#  echo $path_name/$temp
#  echo ${temp%%.*}
  mv $path_name/$temp $path_name/${temp%%.*}.fq.gz
done
#${url##*/} 表示从左边开始删除最后（最右边）一个 / 号及左边的所有字符