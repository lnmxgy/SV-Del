# _*_ coding:utf-8 _*_
#根据输入文件生成执行脚本run1.sh
import getopt
import sys
from collections import OrderedDict

def help():
    print "Usage: python runR1.py ref.fa tumor.sam(bam) nor.sam(bam) k s exlength num /path/picard.jar"

if __name__ == "__main__":
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'k:s:e:n:', ['k=', 's=', 'e=', 'n='])

    imputdict = OrderedDict()
    imputdict['k'] = '10'
    imputdict['s'] = '5'
    imputdict['e'] = '5000'
    imputdict['n'] = '4'
    for i in opts:
        imputdict[i[0][1:]] = i[1]

    nor = args[2].split(".")[0]
    tum = args[1].split(".")[0]
    fp = open("run1.sh", 'w')
    if args[2].split(".")[-1].strip() == "sam":
        fp.write("samtools view -bS " + args[2] + " > " + nor + ".bam &&\n")
        fp.write("samtools  sort  " + nor + ".bam " + nor + ".sort &&\n")
        fp.write("samtools  index  " + nor + ".sort.bam &&\n")
        fp.write("samtools view -bS " + args[1] + " > " + tum + ".bam &&\n")
        fp.write("samtools  sort  " + tum + ".bam " + tum + ".sort &&\n")
        fp.write("samtools  index  " + tum + ".sort.bam &&\n")
    if args[2].split(".")[-1].strip() == "bam":
        fp.write("samtools  sort  " + args[2] + " " + nor + ".sort &&\n")
        fp.write("samtools  index  " + nor + ".sort.bam &&\n")
        fp.write("samtools view -h " + args[2] + " > " + nor + ".sam &&\n")
        fp.write("samtools  sort  " + args[1] + " " + tum + ".sort &&\n")
        fp.write("samtools  index  " + tum + ".sort.bam &&\n")
        fp.write("samtools view -h " + args[1] + " > " + tum + ".sam &&\n")
    fp.write("\ndelly call -t DEL -o s.bcf -g " + args[0] + " " + tum + ".sort.bam " + " " + nor + ".sort.bam &&\n")
    fp.write("svprops s.bcf > s.tab &&\n")
    fp.write("java -Xmx4g -jar " + args[3] + " CollectInsertSizeMetrics I=" + tum + ".bam O=insertsize.txt H=insertsize.pdf &&\n")
    fp.write("python ../main.py " + args[0] + " s.tab " + tum + ".sam " + imputdict['k'] + " " + imputdict['s'] + " ./result1/ ./result2/ " + imputdict['e'] + " " + imputdict['n'])
    fp.close()






