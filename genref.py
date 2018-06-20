# -*- coding: utf-8 -*-
#估计insert size均值和方差。产生第一种情况的ref',断点长度+左右(2u-readlength)bp。并从sam文件中找到含有snv或insert size不对的read

from collections import OrderedDict

def findread(orfile, breakp):
    readlength = 0
    read_seqAll = OrderedDict()
    for breakpi in breakp:
        key = breakpi[0] + "_" + breakpi[1] + ".smvalue"
        read_seqAll[key] = OrderedDict()
    with open(orfile, 'r') as fp:
        while True:
            line1 = fp.readline().strip('\n')
            if not line1:
                break
            if line1[0] == '@':#去掉前面解释的行
                continue
            line2 = fp.readline().strip('\n')
            if not line2:
                break
            tmp1 = line1.split()
            tmp2 = line2.split()

            if tmp1[5] != '150M' or tmp2[5] != '150M':
                i = 0
                while i < len(breakp):#保证取得的read不跨断点
                    if int(tmp1[3]) > int(breakp[i][0]) - 5000 and int(tmp1[3]) < int(breakp[i][0]) + int(breakp[i][1]) + 5000:
                        key = breakp[i][0] + "_" + breakp[i][1] + ".smvalue"
                        read_seqAll[key]['1' + tmp1[0].strip()] = tmp1[9].strip().upper()
                        read_seqAll[key]['2' + tmp2[0].strip()] = tmp2[9].strip().upper()
                        readlength = len(tmp2[9].strip())
                        break
                    if int(tmp2[3]) > int(breakp[i][0]) - 5000 and int(tmp2[3]) < int(breakp[i][0]) + int(breakp[i][1]) + 5000:
                        read_seqAll[key]['1' + tmp1[0].strip()] = tmp1[9].strip().upper()
                        read_seqAll[key]['2' + tmp2[0].strip()] = tmp2[9].strip().upper()
                        break
                    i += 1
    fp.close()
    return read_seqAll, readlength

def readref(reffile):#读原始的100万序列
    totalline = ''
    with open(reffile, 'r') as Fileout:
        Fileout.readline()
        while True:
            line = Fileout.readline().strip('\n')
            if not line:
                break
            totalline += line
    Fileout.close()
    return totalline.upper()

def reddellybrekp(dellyfile):#读delly识别的断点文件
    data = []
    fp = open(dellyfile, 'r')
    fp.readline()
    while True:
        line = fp.readline()
        if not line:
            break
        tmp = line.split()
        datai = [tmp[1], tmp[5]]
        data.append(datai)
    fp.close()
    return data

def guessinsertsize(sam, breakp):
    read_seq, readlength = findread(sam, breakp)
    u= 0
    theta = 0
    with open("insertsize.txt", 'r') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            if line[0 : 6] == "MEDIAN":
                line = fp.readline().strip('\n')
                tmp = line.split()
                u = int(tmp[0])
                theta = int(tmp[-1]) / 2 + 1
    return u, theta, read_seq, readlength

def genref(oref, breakp, insertsize, readlength):#产生新的ref'
    ref = OrderedDict()
    for i in range(0, len(breakp)):
        startp = int(breakp[i][0])
        delen = int(breakp[i][1])
        nref = oref[startp - (2 * insertsize - readlength) : startp + delen + (2 * insertsize - readlength)]
        refkey = breakp[i][0] + '_' + breakp[i][1] + ".smvalue"
        ref[refkey] =  nref
    return ref







