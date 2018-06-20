# _*_ coding:utf-8 _*_
#./main chr19_100w.fa s.tab s150.sam 10 10 ./result1/ ./result2/ ./clustered/
from mapping import *
from breaktext import *
from genref import *
from sh import *
import sys
import datetime
import copy
import numpy

def dealbreakp(breakp, clonf, breakpfl, subCloneReadCount):
    if clonf.has_key(len(breakp)):
        clonf[len(breakp)] = clonf[len(breakp)] + 1
    else:
        clonf[len(breakp)] = 1
    breakpfl.append(breakp)

    tempdict = OrderedDict()
    for key, value in breakp.iteritems():
        tempkey = int(key.split()[-1])
        tempdict[tempkey] = len(value)
    tempdict = OrderedDict(sorted(tempdict.items(), key=lambda e:e[0], reverse=True))#按deletion长度排序，最长的为最后一代克隆
    everysubCloneReadCount = []
    for key, value in tempdict.iteritems():
        everysubCloneReadCount.append(int(value))
    subCloneReadCount.append(everysubCloneReadCount)#依次得到每个断点各个子克隆的deletion的支持度


if __name__ == '__main__':
    print datetime.datetime.now()
    oref = readref(sys.argv[1])
    breakpO = reddellybrekp(sys.argv[2])
    u, threetheta, reads_seq, lengthread = guessinsertsize(sys.argv[3], breakpO)#根据sam文件信息确定insertsize正态分布的均值和方差，以及得到非完美匹配的read

    print "insert size: u = " + str(u) + " 3theta = " + str(threetheta)

    #case1 mapping
    clonf = OrderedDict()#用来存取子克隆数的支持度
    breakfl = []#用于记录所有的deletion信息
    nonemapref = []#用于记录第二种情况的deletion信息
    caseIIbpInfo = OrderedDict()#用来存取预测deletion时得到的第二种情况的断点信息
    subCloneReadCount = []#用于对每一个断点区域各个子克隆read支持数记录
    tolerance = 5
    ref = genref(oref, breakpO, u, lengthread)

    refdel = []
    for k, v in ref.iteritems():
        if len(reads_seq[k]) == 0:
            refdel.append(k)
            continue
        reads_seq_res = OrderedDict()
        seed_table1 = make_seed_table(v, int(sys.argv[4]))
        reads_mapping(reads_seq[k], seed_table1,  int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], k, v, reads_seq_res, 3)

        reads_seq_del = set()
        for kdel in reads_seq[k].keys():#删除
            if kdel[1:] in reads_seq_res.keys() and reads_seq_res[kdel[1:]] == 2:
                reads_seq_del.add(kdel)
        for kdel in reads_seq_del:
            del reads_seq[k][kdel]

    for k in refdel:
        del ref[k]
    for k, v in ref.iteritems():
        breakp, flag = findBreak(sys.argv[6], k, v, u, tolerance)#得到第二种情况的断点信息以及第一种情况的断点
        if flag == 0:
            continue
        halfreflen = (len(v) - (int(k.split(".")[0].split("_")[1]))) / 2
        if len(breakp) < 1:#当断点的个数<=1时为第二中情况
            nonemapref.append(k)
            caseIIbpInfo[k] = "none"
            continue
        elif len(breakp) == 1:
            if breakp.values()[0] > 2:
                if len(breakp.keys()[0].split()) > 2:
                    caseIIbppl = int(breakp.keys()[0].split()[2])
                else:
                    caseIIbppl = int(breakp.keys()[0].split()[1])
                caseIIbpp = int(breakp.keys()[0].split()[0]) + int(k.split("_")[0]) - halfreflen
                caseIIbpInfo[k] = str(caseIIbpp) + " " + str(caseIIbppl)
                nonemapref.append(k)
            else:
                caseIIbpInfo[k] = "none"
            continue
        dealbreakp(breakp, clonf, breakfl, subCloneReadCount)

    #case II mapping
    case2ref = gecase2ref(oref, nonemapref, caseIIbpInfo, u, lengthread, int(sys.argv[8]))
    for k, v in case2ref.iteritems():
        seed_table2 = make_seed_table(v, int(sys.argv[4]))
        reads_mapping(reads_seq[k], seed_table2,  int(sys.argv[4]), int(sys.argv[5]), sys.argv[7], k, v, 0, 0)
    #case II
    for k, v in case2ref.iteritems():
        print k
        snv, readsnv, points, bpp, bppl, halfreflen = Ham(sys.argv[7] + k, k, v, u, threetheta, oref, len(reads_seq[k].values()[0]), reads_seq[k], caseIIbpInfo)
        if halfreflen == 0:
            continue
        count, everysubCloneReadCount = cluster(snv, readsnv, points, lengthread, bpp, bppl, halfreflen, int(sys.argv[9]))
        if count == 0:
            continue
        if clonf.has_key(count):
            clonf[count] = clonf[count] + 1
        else:
            clonf[count] = 1
        subCloneReadCount.append(everysubCloneReadCount)#将第二种情况下子克隆的read信息保留起来


    clonf = OrderedDict(sorted(clonf.items(), key=lambda e:e[1], reverse=True))
    if len(clonf) > 1:
        if clonf.values()[0] - clonf.values()[1] < 10:
            num = max(clonf.keys()[0], clonf.keys()[1])
        else:
            num = clonf.keys()[0]
    else:
        num = clonf.keys()[0]
    print "subclone num: ", num

    print clonf

    selectBreak(nonemapref, copy.deepcopy(ref), breakfl, num)

    total = []
    subCloneReadCountLast = []
    for i in range(0, len(subCloneReadCount)):
        if len(subCloneReadCount[i]) > num:
            subCloneReadCountLast.append(numpy.array(subCloneReadCount[i][:num]))
            total.append(numpy.sum(numpy.array(subCloneReadCount[i][:num])))
            print "*********************"
        elif len(subCloneReadCount[i]) == num:
            subCloneReadCountLast.append(numpy.array(subCloneReadCount[i]))
            total.append(numpy.sum(numpy.array(subCloneReadCount[i])))

    for i in range(0, len(subCloneReadCountLast)):
        subCloneReadCountLast[i] = subCloneReadCountLast[i] * 1.0 / total[i]

    ratio = map(sum, zip(*subCloneReadCountLast))
    ratio = numpy.array(ratio)
    print ratio / len(subCloneReadCountLast)

    print datetime.datetime.now()







