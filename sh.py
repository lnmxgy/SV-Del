# _*_ coding:utf-8 _*_
#二次哈希，将150长的read分为15个单位，每个求一个相似度的值，将所有相似度的值求平均，使用read文件中相对位置信息确定read的比对位置，不进行mapping
#check
from mapping import *
from collections import defaultdict
import copy
from breaktext import *
import math

def find(shresult, listBoth, listOne):
    with open(shresult, 'r') as fp:
        fp.readline()
        while True:
            line1 = fp.readline().strip('\n')
            line2 = fp.readline().strip('\n')
            if not line1:
                break
            if not line2:
                break
            tmp1 = line1.split("|")
            tmp2 = line2.split("|")
            if len(tmp1) > 4 and len(tmp2) > 4:
                temp1 = tmp1[3].split("_")
                temp2 = tmp2[3].split("_")
                if int(temp1[1]) < int(temp2[1]):#如果两端都有比对结果，将汉明距离大的一端写入
                    listi = [tmp2[0].strip(), tmp2[1].strip(), tmp2[2].strip(), tmp2[3].strip()]
                    listBoth.append(listi)
                else:
                    listi = [tmp1[0].strip(), tmp1[1].strip(), tmp1[2].strip(), tmp1[3].strip()]
                    listBoth.append(listi)
            else:
                if len(tmp1) > 4:#如果只有一端有比对结果，将没有比对结果的一端及有比对结果的一端比对位置写入文件
                    listi = [tmp2[0].strip(), tmp2[1].strip(), tmp2[2].strip(), tmp1[3].strip()]
                    listOne.append(listi)
                if len(tmp2) > 4:
                    listi = [tmp1[0].strip(), tmp1[1].strip(), tmp1[2].strip(), tmp2[3].strip()]
                    listOne.append(listi)
    fp.close()

def Ham(shresult, refk, refseq, u, theta, oref, readlength, reads_seq, caseIIbpInfo):
    bpp = int(refk.split(".")[0].split("_")[0])#断点位置
    bppl = int(refk.split(".")[0].split("_")[1])#断点长度
    listBoth = []#用于存储两端都有比对信息汉明距离大的一端的read信息
    listOne = []#用于存储只有一端有比对信息中没有比对信息的一端read
    snv = defaultdict(list)#用于保存每个snv包含的read信息
    readsnv = defaultdict(list)#用于保存每个read包含snv的信息
    points = OrderedDict()#用于保留read相对于refmn的比对位置
    halfreflen = len(refseq) / 2

    find(shresult, listBoth, listOne)#读与该ref匹配上的read文件
    if len(listBoth) == 0 and len(listOne) == 0:
        print "mapping exist empty file", refk
        return snv, readsnv, points, 0, 0, 0

    readPos = dict()
    if caseIIbpInfo[refk] == "none":
        Tmn = dict()
        for i in range(0, len(listBoth)):#使用listBoth中的read推断出真正的断点
            canloc = int(listBoth[i][3].split("_")[0])
            if canloc < halfreflen - readlength or canloc > halfreflen:#只考虑能跨越断点的read
                continue
            flag = listBoth[i][1]
            if 'z' in flag:
                read = listBoth[i][2]
            else:
                read = ''.join(["ATCG"["TAGC".index(x)] for x in listBoth[i][2].strip()[::-1]])
            listBoth[i][2] = read
            reads_seq[listBoth[i][0]] = read
            mn = dict()
            for m in range(bpp - 10, bpp + 11):#改变delly预测的断点构建新的断点
                for n in range(bppl - 10, bppl + 11):
                    if m + n + halfreflen > len(oref):#越界
                        return snv, readsnv, points, 0, 0, 0
                    refmn = oref[m - 200 : m] + oref[m + n : m + n + 200]
                    can = findLongestCommonSubtringWithTolerance(read, refmn, 2)
                    mdis = readlength
                    frefstart = 0
                    freadstart = 0
                    for cani in can:
                        mmax, readstart, refstart = cani
                        dis = seed_distance(read, refmn[refstart - readstart : refstart - readstart + readlength])
                        if dis < mdis:
                            mdis = dis
                            frefstart = refstart
                            freadstart = readstart
                    keymn = str(m) + "_" + str(n) + " " + str(listBoth[i][0].strip())
                    mn[keymn] = mdis#记录每一条read和每一条参考序列的最小值
                    readPos[keymn] = frefstart - freadstart#记录每一条read和每一条参考序列的最小值时read在该参考序列上的位置
            mn = OrderedDict(sorted(mn.items(), key=lambda e:e[1], reverse=False))#每一条read与多条参考序列的最小值
            mindis = mn.values()[0]
            for mapk, mapv in mn.iteritems():#每一条参考序列
                if mapv == mindis:
                    keyT = mapk.split()[0]
                    if Tmn.has_key(keyT):
                        Tmn[keyT] += 1
                    else:
                        Tmn[keyT] = 1

        Tmn = OrderedDict(sorted(Tmn.items(), key=lambda e:e[1], reverse=True))
        if len(Tmn) == 0:
            print "no appropriate read"
            return snv, readsnv, points, 0, 0, 0
        m = int(Tmn.keys()[0].split("_")[0])#选出支持度最高的ref'
        n = int(Tmn.keys()[0].split("_")[1])
    else:
        m = int(caseIIbpInfo[refk].split()[0])
        n = int(caseIIbpInfo[refk].split()[1])
    print "deltion : ", m, n
    refmn = oref[m - halfreflen : m] + oref[m + n : m + n + halfreflen]

    for i in range(0, len(listBoth)):
        flag = listBoth[i][1]
        if 'z' in flag:
            read = listBoth[i][2]
        else:
            read = ''.join(["ATCG"["TAGC".index(x)] for x in listBoth[i][2].strip()[::-1]])
        listBoth[i][2] = read
        reads_seq[listBoth[i][0]] = read
        canloc = int(listBoth[i][3].split("_")[0])
        if canloc < 500 or canloc > 9500:#只考虑能跨越断点的read
            continue
        if caseIIbpInfo[refk] == "none":
            keyread = str(m) + "_" + str(n) + " " + str(listBoth[i][0].strip())
            if keyread in readPos.keys():
                minpos = readPos[keyread] + halfreflen - 200
            else:
                continue
        else:
            minpos = canloc
        points[listBoth[i][0]] = minpos
        for posi in range(0, readlength):
            if read[posi] != refmn[minpos]:
                snv[minpos].append(listBoth[i][0])
                readsnv[listBoth[i][0]].append(minpos)
            minpos += 1

    for i in range(0, len(listOne)):
        othermloc = int(listOne[i][3].split("_")[0])
        if othermloc < 500 or othermloc > 9500:
            continue
        flag = listOne[i][1]
        mindis = 150
        minpos = 0
        readf = ''
        if othermloc > halfreflen:
            canloc = othermloc - (u - readlength)
        else:
            canloc = othermloc + (u - readlength)
        if int(flag[1 : ]) > 3:
            if 'z' in flag:
                read = listOne[i][2]
            else:
                read = ''.join(["ATCG"["TAGC".index(x)] for x in listOne[i][2].strip()[::-1]])
            for start in range(canloc - theta, canloc + theta):
                dis = seed_distance(read, refmn[start : start + readlength])
                if(dis < mindis):
                    mindis = dis
                    minpos = start
                    readf = read
        else:
            readz = listOne[i][2]
            readv = ''.join(["ATCG"["TAGC".index(x)] for x in listOne[i][2].strip()[::-1]])
            for start in range(canloc - theta, canloc + theta):
                disz = seed_distance(readz, refmn[start : start + readlength])
                disv = seed_distance(readv, refmn[start : start + readlength])
                if(disz < disv):
                    read = readz
                else:
                    read = readv
                if(min(disz, disv) < mindis):
                    mindis = min(disz, disv)
                    minpos = start
                    readf = read
        listOne[i][2] = readf
        reads_seq[listOne[i][0]] = readf
        if mindis > 60:
            continue
        points[listOne[i][0]] = minpos
        for posi in range(0, readlength):
            if readf[posi] != refmn[minpos]:
                snv[minpos].append(listOne[i][0])
                readsnv[listOne[i][0]].append(minpos)
            minpos += 1


    readdel = set()#存测序错误的read

    for k in snv.keys():
        snvdict = defaultdict(list)#存每一个变异位点的变异情况,用包含该位点的read对该位点进行投票
        cov = 0
        for x in readsnv.keys():
            if k >= points[x] and k < points[x] + readlength:
                cov += 1
                snvdict[reads_seq[x][k - points[x]]].append(x)
        snvdict = OrderedDict(sorted(snvdict.items(), key=lambda e:len(e[1]), reverse=True))
        if list(snvdict.keys())[0] != refmn[k] and len(list(snvdict.values())[0]) > 1:
            for point in snvdict.keys():
                if point not in [list(snvdict.keys())[0], refmn[k]]:
                    for rn in snvdict[point]:
                        readdel.add(rn)

        elif len(snvdict.keys()) >= 2 and len(list(snvdict.values())[1]) > 1:#第二种碱基支持度大于5
            for point in snvdict.keys():
                if point not in [list(snvdict.keys())[1], refmn[k]]:#删除测序错误的read
                    for rn in snvdict[point]:
                        readdel.add(rn)

    for readtodel in readdel:#删除read时保持snv和readsnv一致
        del readsnv[readtodel]
    for k in snv.keys():
        snv[k] = list(set(snv[k]).difference(set(readdel)))
        snv[k].sort()

    snvSave = []
    for mm, nn in snv.iteritems():
        if len(nn) > 3:#
            snvSave.append(mm)

    if len(snvSave) < 5:
        print "simple deletion", m, n
        return snv, readsnv, points, 0, 0, 0


    snvamm = list(snv.keys())
    snvamm.sort()
    with open("snvall.txt", 'a') as fl:
        for vc in snvamm:
            if len(snv[vc]) > 2:
                if vc > halfreflen:
                    fl.write(str(vc + m - halfreflen + n + 1) + " ")
                    fl.write(str(len(snv[vc])) + '\n')
                else:
                    fl.write(str(vc + m - halfreflen + 1) + " ")
                    fl.write(str(len(snv[vc])) + '\n')
    fl.close()



    snvSave.sort()
    print snvSave

    left = snvSave[0]
    right = snvSave[len(snvSave) - 1]
    lenhighdensity = right - left
    print lenhighdensity
    mutationrate = len(snvSave) * 1.0 / lenhighdensity
    print mutationrate
    if mutationrate > 0.05:
        while mutationrate > 0.05:
            lenhighdensity += 10
            mutationrate = len(snvSave) * 1.0 / lenhighdensity
        print mutationrate, lenhighdensity
    else:
        while mutationrate < 0.05:
            lenhighdensity -= 10
            mutationrate = len(snvSave) * 1.0 / lenhighdensity
        print mutationrate, lenhighdensity

    snvdel = list(set(snv.keys()).difference(set(snvSave)))
    readdel = set()
    snv, readsnv, snvdel, readdel = delsnv(snv, readsnv, snvdel, readdel, 0)
    return snv, readsnv, points, m, n, halfreflen



def delsnv(snv, readsnv, snvdel, readdel, snvtolerance):
    for k in snv.keys():#去掉支持度小的snv保持snv和readsnv一致
        if len(snv[k]) < snvtolerance:
            snvdel.append(k)
    for k in snv.keys():
        if k in snvdel:
            del snv[k]
    for k in readsnv.keys():#从read中删除上一步去掉的snv
        readsnv[k] = list(set(readsnv[k]).difference(set(snvdel)))
        readsnv[k].sort()

    for k in readsnv.keys():
        if len(readsnv[k]) == 0:
            readdel.add(k)
    for readtodel in readdel:#删除read时保持snv和readsnv一致
        del readsnv[readtodel]
    for k in snv.keys():
        snv[k] = list(set(snv[k]).difference(set(readdel)))
        snv[k].sort()

    return snv, readsnv, snvdel, readdel


def cluster(snvo, readsnv, points, readlength, bpp, bppl, halfreflen, num):
    if len(snvo) == 0:
        print "snv 0"
    if len(readsnv) == 0:
        print "readsnv 0"
    freadsnv = copy.deepcopy(readsnv)
    snv = copy.deepcopy(snvo)
    snveveryclone = []
    readeveryclone = []
    count = 0
    snvall = list(snv.keys())
    snvall.sort()
    while(len(snvall) > 0):
        print "**********************"
        if count > 10:
            break
        count += 1
        snv4read = []

        hasnosnvdic = defaultdict(set)
        for key, valu in freadsnv.iteritems():
            minpos = 0
            maxpos = 0
            for snvpos in range(0, len(snvall)):
                if points[key] <= snvall[snvpos]:
                    minpos = snvpos
                    break
            for snvpos in range(0, len(snvall)):
                if points[key] + readlength - 1 < snvall[snvpos]:
                    maxpos = snvpos
                    break
            if minpos != 0 and maxpos == 0:
                maxpos = len(snvall)
            shouldhassnv = snvall[minpos : maxpos]
            hasnosnv = [v for v in shouldhassnv if v not in valu]
            hasnosnvdic[key] = set(hasnosnv)
            if len(hasnosnv) == 0:
                snv4read.append(key)

        readdel = set(snv4read)
        snvdel = []
        snv, freadsnv, snvdel, readdel = delsnv(snv, freadsnv, snvdel, readdel, 0)#删read
        readdel = set()
        snv, freadsnv, snvdel, readdel = delsnv(snv, freadsnv, snvdel, readdel, num)#删snv


        snv4read = []
        for key, valu in freadsnv.iteritems():
            minpos = 0
            maxpos = 0
            for snvpos in range(0, len(snvall)):
                if points[key] <= snvall[snvpos]:
                    minpos = snvpos
                    break
            for snvpos in range(0, len(snvall)):
                if points[key] + readlength - 1 < snvall[snvpos]:
                    maxpos = snvpos
                    break
            if minpos != 0 and maxpos == 0:
                maxpos = len(snvall)
            shouldhassnv = snvall[minpos : maxpos]
            hasnosnv = [v for v in shouldhassnv if v not in valu]
            hasnosnvdic[key] = set(hasnosnv)
            if len(hasnosnv) == 0:
                snv4read.append(key)

        if len(snv4read) < num:
            readdel = set(snv4read)
            snv, freadsnv, snvdel, readdel = delsnv(snv, freadsnv, snvdel, readdel, 0)
            snv, freadsnv, snvdel, readdel = delsnv(snv, freadsnv, snvdel, readdel, 1)#删snv

        snvall = list(snv.keys())
        snvall = list(set(snvall).difference(set(snvdel)))
        snvall.sort()#返回值为空
        snvdel.sort()
        snveveryclone.append(snvdel)
        readeveryclone.append(len(snv4read))
        print snvall
        print snvdel
        print len(snv4read)

    if count >= 10:
        return 0, 0

    with open('snv.txt', 'a') as fp:
        fp.write("*" + str(bpp) + " " + str(bppl) + '\n')
        for m in range(0, len(snveveryclone)):
            fp.write("subclone " + str(m + 1) + ":")
            for n in snveveryclone[len(snveveryclone) - m - 1]:
                if n > halfreflen:
                    fp.write(str(n + bpp - halfreflen + bppl + 1) + ' ')
                else:
                    fp.write(str(n + bpp - halfreflen + 1) + ' ')
            fp.write('\n')
        fp.write("$\n")
    return len(snveveryclone), readeveryclone


